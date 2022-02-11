#include <string>
#include <tuple>

#include "../lib/commandIO/src/commandIO.h"
#include "../lib/fastp/src/writer.h"
#include "../lib/trie/src/trie.tcc"

#include "cluster.h"
#include "fastq.h"
#include "leaf.h"
#include "log.h"

using std::ios;

/*!
 */
tuple<size_t, size_t> readData(
    Trie<4, NLeaf>& trie, vector<string> files, size_t length,
    ofstream& log) {
  time_t start = startMessage(log, "Reading data");
  size_t total = 0;
  size_t line = 0;
  for (vector<Read*> reads: readFiles(files)) {
    Word word = makeWord(reads, length);
    if (!word.filtered) {
      NLeaf* leaf = trie.add(word.data);
      leaf->lines.push_back(line++);
    }
    total++;
  }
  endMessage(log, start);

  return tuple<size_t, size_t>(total, line);
}

/*!
 */
size_t findNeighbours(Trie<4, NLeaf>& trie, size_t distance, ofstream& log) {
  size_t start = startMessage(log, "Calculating neighbours");
  size_t nonDuplicates = 0;
  for (Result<NLeaf> walkResult: trie.walk()) {
    for (Result<NLeaf> hammingResult: trie.hamming(
        walkResult.path, distance)) {
      if (walkResult.leaf != hammingResult.leaf) {
        walkResult.leaf->neighbours.push_back(hammingResult.leaf);
        hammingResult.leaf->neighbours.push_back(walkResult.leaf);
      }
    }
    nonDuplicates++;
  }
  endMessage(log, start);

  return nonDuplicates;
}

/*!
 */
vector<Cluster*> findClusters(Trie<4, NLeaf>& trie, ofstream& log) {
  size_t start = startMessage(log, "Calculating clusters");
  vector<Cluster*> clusters;
  size_t id = 0;
  for (Result<NLeaf> result: trie.walk()) {
    if (!result.leaf->cluster) {
      Cluster* cluster = new Cluster(id++);
      assignCluster(result.leaf, cluster);
      clusters.push_back(cluster);
    }
  }
  endMessage(log, start);

  return clusters;
}

/*!
 */
void writeResults(
    Trie<4, NLeaf>& trie, vector<string> files, size_t length,
    string dirName, ofstream& log) {
  size_t start = startMessage(log, "Writing results");
  vector<Writer*> outFiles;
  Options options;
  for (string name: makeFileNames(files, dirName)) {
    outFiles.push_back(new Writer(&options, name, options.compression));
  }
  for (vector<Read*> reads: readFiles(files)) {
    Word word = makeWord(reads, length);
    if (!word.filtered) {
      Node<4, NLeaf>* node = trie.find(word.data);
      if (
          !node->leaf->cluster->visited &&
          node->leaf->cluster->maxLeaf == node->leaf) {
        for (size_t i = 0; i < reads.size(); i++) {
          string s = reads[i]->toString();
          outFiles[i]->write(s.c_str(), s.size());
        }
        node->leaf->cluster->visited = true;
      }
    }
  }
  for (Writer* w: outFiles) {
    delete w;
  }
  endMessage(log, start);
}

/*!
 */
tuple<map<size_t, size_t>, map<size_t, size_t>> runStatistics(
    Trie<4, NLeaf>& trie, ofstream& log) {
  size_t start = startMessage(log, "Calculating count and neighbour stats");

  map<size_t, size_t> counts;
  map<size_t, size_t> neighbours;
  for (Result<NLeaf> result: trie.walk()) {
    counts[result.leaf->count]++;
    neighbours[result.leaf->neighbours.size()]++;
  }
  endMessage(log, start);

  return tuple<map<size_t, size_t>, map<size_t, size_t>>(
    counts, neighbours);
}

/*!
 */
void writeStatistics(
    map<size_t, size_t> counts, map<size_t, size_t> neighbours,
    map<size_t, size_t> clusters, size_t total, size_t useable, size_t unique,
    size_t clusterSize, string dirName) {
  ofstream output(addDir("counts.dat", dirName), ios::out | ios::binary);
  for (pair<size_t, size_t> count: counts) {
    output << count.first << ' ' << count.second << '\n';
  }
  output.close();

  output.open(addDir("neigh.dat", dirName), ios::out | ios::binary);
  for (pair<size_t, size_t> count: neighbours) {
    output << count.first << ' ' << count.second << '\n';
  }
  output.close();

  output.open(addDir("clusters.dat", dirName), ios::out | ios::binary);
  for (pair<size_t, size_t> count: clusters) {
    output << count.first << ' ' << count.second << '\n';
  }
  output.close();

  output.open(addDir("stats.dat", dirName), ios::out | ios::binary);
  output << "total: " << total << '\n';
  output << "useable: " << useable << '\n';
  output << "unique: " << unique << '\n';
  output << "clusters: " << clusterSize << '\n';
  output.close();
}

/*!
 * Determine duplicates.
 *
 * \param wordLength Read length.
 * \param distance Maximum hamming distance between reads.
 * \param logName Log file.
 * \param dirName Output directory.
 * \param files FastQ files.
 */
void dedup(
    size_t wordLength, size_t distance, string logName, string dirName,
    bool runStats, bool write, vector<string> files) {
  Trie<4, NLeaf> trie;
  size_t length = wordLength / files.size();

  ofstream log(logName.c_str(), ios::out | ios::binary);

  tuple<size_t, size_t> input = readData(trie, files, length, log);
  size_t nonDuplicates = findNeighbours(trie, distance, log);
  vector<Cluster*> clusters = findClusters(trie, log);

  if (write) {
    writeResults(trie, files, length, dirName, log);
  }
  if (runStats) {
    tuple<map<size_t, size_t>, map<size_t, size_t>> stats = runStatistics(
      trie, log);
    map<size_t, size_t> cStats = clusterStats(clusters);
    writeStatistics(
      get<0>(stats), get<1>(stats), cStats, get<0>(input), get<1>(input),
      nonDuplicates, clusters.size(), dirName);
  }

  log.close();

  freeClusters(clusters);
}


/*
 * Argument parsing.
 */
int main(int argc, char* argv[]) {
  CliIO io(argc, argv);

  interface(
    io,
    dedup, argv[0], "Deduplicate a dataset.", 
      param("-l", 24, "word length"),
      param("-d", 1, "distance"),
      param("-o", "/dev/stderr", "log file name"),
      param("-f", ".", "output directory"),
      param("-s", false, "run statistics"),
      param("-w", true, "do not write output"),
      param("files", "FastQ files"));

  return 0;
}
