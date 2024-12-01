#include <filesystem>
#include <tuple>

#include "../lib/commandIO/src/commandIO.h"
#include "../lib/fastp/src/writer.h"
#include "../lib/trie/src/trie.tcc"

#include "cluster.h"
#include "fastq.h"
#include "leaf.h"
#include "log.h"

using std::filesystem::create_directories;
using std::ios;
using std::tie;

/*! Peek at the header of the first read, to determine the size of the UMI, if
 * any.
 *
 * \param filename Input file name.
 *
 * \return Size of the UMI in the header.
 */
size_t peekUMI(string const filename) {
  FastqReader reader {filename.c_str()};
  Read* read {reader.read()};

  size_t umiSize {extractUMI(read).size()};

  delete read;

  return umiSize;
}

/*! Pre-compute the nucleotides to take from the UMI header, and from each of
 * the input files.
 */
tuple<size_t, vector<size_t>> preCompute(
    vector<string> const files, size_t const wordLength) {
  // Peek at the header of the first read in the first file to get the UMI size.
  size_t headerUMISize {peekUMI(files.front())};

  // Ensure we do not take a negative amount of nucleotides from the files.
  size_t fromFile {0};
  if (wordLength > headerUMISize) {
    fromFile = wordLength - headerUMISize;
  }

  // Calculate how many nucleotides to take from each read. Any remainder will
  // be taken from the last file.
  vector<size_t> ntToTake {ntFromFile(files.size(), fromFile)};

  // Ensure we do not take more than `wordLength` from the UMI header.
  if (wordLength < headerUMISize) {
    headerUMISize = wordLength;
  }

  return tuple<size_t, vector<size_t>>(headerUMISize, ntToTake);
}

/*! Populate a trie with words extracted from FastQ files.
 *
 * \param trie Trie.
 * \param files Input file names.
 * \param wordLength Word length.
 * \param log Log handle.
 *
 * \return Total and usable number of reads.
 */
tuple<size_t, size_t> readData(
    Trie<4, NLeaf>& trie, vector<string> const files, size_t const wordLength,
    ofstream& log) {

  // Pre calculate some values so that we do not have to re-calculate them for
  // every single read.
  size_t headerUMISize;
  vector<size_t> ntToTake;
  tie(headerUMISize, ntToTake) = preCompute(files, wordLength);

  time_t nt_start {startMessage(log, "Determing nucleotides to take")};
  endMessage(log, nt_start);

  log << "  header: " << headerUMISize;
  for (size_t i {0}; i < ntToTake.size(); ++i) {
    log << "\n  " << files[i] << ": " << ntToTake[i];
  }
  log << "\n";

  time_t start {startMessage(log, "Reading data")};
  size_t total {0};
  size_t usable {0};
  for (vector<Read*> const& reads: readFiles(files)) {
    Word word {makeWord(reads, ntToTake, headerUMISize)};
    if (not word.filtered) {
      trie.add(word.data);
      usable++;
    }
    total++;
  }
  endMessage(log, start);

  return tuple<size_t, size_t>(total, usable);
}

/*! Calculate neighbours for every word in a trie.
 *
 * \param trie Trie.
 * \param distance Maximum neighbour distance.
 * \param log Log handle.
 *
 * \return Number of unique words.
 */
size_t findHammingNeighbours(
    Trie<4, NLeaf> const& trie, size_t const distance, ofstream& log) {
  time_t start {startMessage(log, "Calculating neighbours using Hamming distance")};
  size_t unique {0};
  for (Result<NLeaf> const& walkResult: trie.walk()) {
    for (Result<NLeaf> const& hammingResult: trie.asymmetricHamming(
        walkResult.path, distance)) {
      if (walkResult.leaf != hammingResult.leaf) {
        walkResult.leaf->neighbours.push_back(hammingResult.leaf);
        hammingResult.leaf->neighbours.push_back(walkResult.leaf);
      }
    }
    unique++;
  }
  endMessage(log, start);

  return unique;
}

/*! Calculate neighbours for every word in a trie.
 *
 * \param trie Trie.
 * \param distance Maximum neighbour distance.
 * \param log Log handle.
 *
 * \return Number of unique words.
 */
size_t findEditNeighbours(
    Trie<4, NLeaf> const& trie, size_t const distance, ofstream& log) {
  time_t start {startMessage(log, "Calculating neighbours using Levenshtein distance")};
  size_t unique {0};

  for (Result<NLeaf> const& walkResult: trie.walk()) {
    for (Result<NLeaf> const& editResult: trie.asymmetricLevenshtein(
        walkResult.path, distance)) {
      if (walkResult.leaf != editResult.leaf) {
        walkResult.leaf->neighbours.push_back(editResult.leaf);
        editResult.leaf->neighbours.push_back(walkResult.leaf);
      }
    }
    unique++;
  }
  endMessage(log, start);

  return unique;
}

/*! Group neighbours into clusters.
 *
 * \param trie Trie.
 * \param log Log handle.
 *
 * \return Number of clusters.
 */
vector<Cluster*> findClusters(
    Trie<4, NLeaf>& trie, bool const maximum, ofstream& log) {
  time_t start{};
  if (maximum) {
    start = startMessage(log, "Calculating maximum clusters");
  }
  else {
    start = startMessage(log, "Calculating directional clusters");
  }
  vector<Cluster*> clusters;
  size_t id {0};
  for (Result<NLeaf> const& result: trie.walk()) {
    if (not result.leaf->cluster) {
      Cluster* cluster {new Cluster {id++}};
      if (maximum) {
        assignMaxCluster(result.leaf, cluster);
      }
      else {
        assignDirectionalCluster(result.leaf, cluster);
      }
      clusters.push_back(cluster);
    }
  }
  endMessage(log, start);

  return clusters;
}

/*! Filter FastQ files for duplicates.
 *
 * \param trie Trie.
 * \param files Input file names.
 * \param wordLength Word length.
 * \param dirName Output directory.
 * \param log Log handle.
 */
void writeFiltered(
    Trie<4, NLeaf> const& trie, vector<string> const files,
    size_t const wordLength, string const dirName, ofstream& log) {
  time_t start {startMessage(log, "Writing filtered results")};

  // Pre calculate some values so that we do not have to re-calculate them for
  // every single read.
  size_t headerUMISize;
  vector<size_t> ntToTake;
  tie(headerUMISize, ntToTake) = preCompute(files, wordLength);

  vector<Writer*> outFiles;
  Options options;
  for (string const& name: makeFileNames(files, dirName, "dedup")) {
    outFiles.push_back(new Writer(&options, name, options.compression));
  }

  for (vector<Read*> const& reads: readFiles(files)) {
    Word word {makeWord(reads, ntToTake, headerUMISize)};
    if (not word.filtered) {
      Node<4, NLeaf>* node {trie.find(word.data)};
      if (
          !node->leaf->cluster->visited &&
          node->leaf->cluster->maxLeaf == node->leaf) {
        for (size_t i {0}; i < reads.size(); i++) {
          string s {reads[i]->toString()};
          outFiles[i]->write(s.c_str(), s.size());
        }
        node->leaf->cluster->visited = true;
      }
    }
  }

  for (Writer* const w: outFiles) {
    delete w;
  }

  endMessage(log, start);
}

/*! Annotate FastQ files with cluster IDs.
 *
 * \param trie Trie.
 * \param files Input file names.
 * \param wordLength Word length.
 * \param dirName Output directory.
 * \param log Log handle.
 */
void writeAnnotated(
    Trie<4, NLeaf> const& trie, vector<string> const files, size_t const
    wordLength, string const dirName, ofstream& log) {
  time_t start {startMessage(log, "Writing annotated results")};

  // Pre calculate some values so that we do not have to re-calculate them for
  // every single read.
  size_t headerUMISize;
  vector<size_t> ntToTake;
  tie(headerUMISize, ntToTake) = preCompute(files, wordLength);

  vector<Writer*> outFiles;
  Options options;
  for (string const& name: makeFileNames(files, dirName, "annotated")) {
    outFiles.push_back(new Writer(&options, name, options.compression));
  }

  for (vector<Read*> const& reads: readFiles(files)) {
    Word word {makeWord(reads, ntToTake, headerUMISize)};
    if (not word.filtered) {
      Node<4, NLeaf>* node {trie.find(word.data)};

      for (size_t i {0}; i < reads.size(); i++) {
        *reads[i]->mName += ':' + to_string(node->leaf->cluster->id);
        string s {reads[i]->toString()};
        outFiles[i]->write(s.c_str(), s.size());
      }
    }
  }

  for (Writer* const w: outFiles) {
    delete w;
  }

  endMessage(log, start);
}

/*! Make histograms of the number of perfect and nonperfect duplicates.
 *
 * \param trie Trie.
 * \param log Log handle.
 *
 * \return Duplicate statistics histograms.
 */
tuple<map<size_t, size_t>, map<size_t, size_t>> runStatistics(
    Trie<4, NLeaf> const& trie, ofstream& log) {
  time_t start {startMessage(log, "Calculating count and neighbour stats")};

  map<size_t, size_t> counts;
  map<size_t, size_t> neighbours;
  for (Result<NLeaf> const& result: trie.walk()) {
    counts[result.leaf->count]++;
    neighbours[result.leaf->neighbours.size()]++;
  }
  endMessage(log, start);

  return tuple<map<size_t, size_t>, map<size_t, size_t>>(
    counts, neighbours);
}

/*! Write statistics to files.
 *
 * \param counts Perfect duplicate histogram.
 * \param neighbours Nonperfect duplicate histogram.
 * \param clusters 
 * \param total 
 * \param usable 
 * \param unique 
 * \param clusterSize 
 * \param dirName Output directory.
 */
void writeStatistics(
    map<size_t, size_t> const& counts, map<size_t, size_t> const& neighbours,
    map<size_t, size_t> const& clusters, size_t const total,
    size_t const usable, size_t const unique, size_t const clusterSize,
    string const dirName) {
  ofstream output(addDir("counts.dat", dirName), ios::out | ios::binary);
  for (pair<size_t const, size_t> const& count: counts) {
    output << count.first << ' ' << count.second << '\n';
  }
  output.close();

  output.open(addDir("neigh.dat", dirName), ios::out | ios::binary);
  for (pair<size_t const, size_t> const& count: neighbours) {
    output << count.first << ' ' << count.second << '\n';
  }
  output.close();

  output.open(addDir("clusters.dat", dirName), ios::out | ios::binary);
  for (pair<size_t const, size_t> const& count: clusters) {
    output << count.first << ' ' << count.second << '\n';
  }
  output.close();

  output.open(addDir("stats.dat", dirName), ios::out | ios::binary);
  output << "total: " << total << '\n';
  output << "usable: " << usable << '\n';
  output << "unique: " << unique << '\n';
  output << "clusters: " << clusterSize << '\n';
  output.close();
}

/*! Determine duplicates.
 *
 * \param wordLength Read length.
 * \param distance Maximum distance between reads.
 * \param logName Log file.
 * \param dirName Output directory.
 * \param runStats
 * \param write
 * \param files FastQ files.
 */
void humid(
    size_t const wordLength, size_t const distance, string const logName,
    string const dirName, bool const runStats, bool const filter,
    bool const annotate, bool const edit, bool const maximum,
    vector<string> const files) {
  Trie<4, NLeaf> trie;

  ofstream log(logName.c_str(), ios::out | ios::binary);

  tuple<size_t, size_t> input {readData(trie, files, wordLength, log)};

  size_t unique;
  if (edit) {
    unique = findEditNeighbours(trie, distance, log);
  }
  else {
    unique = findHammingNeighbours(trie, distance, log);
  }

  vector<Cluster*> clusters {findClusters(trie, maximum, log)};

  create_directories(dirName);
  if (filter) {
    writeFiltered(trie, files, wordLength, dirName, log);
  }
  if (annotate) {
    writeAnnotated(trie, files, wordLength, dirName, log);
  }
  if (runStats) {
    tuple<map<size_t, size_t>, map<size_t, size_t>> stats {runStatistics(
      trie, log)};
    map<size_t, size_t> cStats {clusterStats(clusters)};
    writeStatistics(
      get<0>(stats), get<1>(stats), cStats, get<0>(input), get<1>(input),
      unique, clusters.size(), dirName);
  }

  log.close();

  freeClusters(clusters);
}


/* Argument parsing. */
int main(int argc, char* argv[]) {
  CliIO io(argc, argv);

  interface(
    io,
    humid, argv[0], "Deduplicate a dataset.",
      param("-n", 24, "word length"),
      param("-m", 1, "allowed mismatches"),
      param("-l", "/dev/stderr", "log file name"),
      param("-d", ".", "output directory"),
      param("-s", false, "calculate statistics"),
      param("-q", true, "write deduplicated FastQ files"),
      param("-a", false, "write annotated FastQ files"),
      param("-e", false, "use edit distance"),
      param("-x", false, "use maximum clustering method"),
      param("files", "FastQ files"));
}
