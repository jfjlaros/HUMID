#include <string>

#include "../lib/commandIO/src/commandIO.h"

#include "../../src/trie.tcc"
#include "../lib/ngs/ngs.h"

using std::ios;
using std::ofstream;
using std::string;

/*!
 * Cluster structure.
 */
struct Cluster {
  size_t id;
  bool visited = false;

  Cluster(size_t id) { this->id = id; }
};

/*!
 * Leaf structure for neighbour finding.
 */
struct NLeaf : Leaf {
  vector<size_t> lines;
  vector<NLeaf*> neighbours;
  Cluster* cluster = NULL;
};


/*!
 * Traverse neighbours to assign cluster IDs.
 *
 * \param leaf Leaf node.
 * \param cluster Cluster.
 */
void assignCluster(NLeaf* leaf, Cluster* cluster) {
  leaf->cluster = cluster;
  for (NLeaf* neighbour: leaf->neighbours) {
    if (!neighbour->cluster) {
      assignCluster(neighbour, cluster);
    }
  }
}

/*!
 * Write a task start message to a log.
 *
 * \param log Log file.
 * \param message Message.
 *
 * \return Task start time.
 */
time_t startMessage(ofstream& log, char const message[]) {
  log << message << "... ";
  log.flush();

  return time(NULL);
}

/*!
 * Write a task end message to a log.
 *
 * \param log Log file.
 * \param start Task start time.
 */
void endMessage(ofstream& log, time_t start) {
  unsigned int seconds = (unsigned int)difftime(time(NULL), start);
  log << "done. (" << seconds / 60 << 'm' << seconds % 60 << "s)\n";
  log.flush();
}

/*!
 * Determine duplicates.
 *
 * \param read1 FastQ file for read 1.
 * \param read2 FastQ file for read 2.
 * \param umi FastQ file for the UMI.
 * \param length Read length.
 * \param distance Maximum hamming distance between reads.
 * \param logName Log file.
 */
void dedup(
    string read1, string read2, string umi, size_t length,
    string r1, string r2, string u,
    size_t distance, string logName) {
  vector<string> files = {read1, read2, umi};
  Trie<4, NLeaf> trie;

  ofstream log(logName.c_str(), ios::out | ios::binary);
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

  start = startMessage(log, "Calculating neighbours");
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

  start = startMessage(log, "Calculating clusters");
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

  start = startMessage(log, "Writing results");
  ofstream r1Out(r1.c_str(), ios::out | ios::binary);
  ofstream r2Out(r2.c_str(), ios::out | ios::binary);
  ofstream uOut(u.c_str(), ios::out | ios::binary);
  for (vector<Read*> reads: readFiles(files)) {
    Word word = makeWord(reads, length);
    if (!word.filtered) {
      Node<4, NLeaf>* node = trie.find(word.data);
      if (!node->leaf->cluster->visited) {
        // Emit reads.
        reads[0]->printFile(r1Out);
        reads[1]->printFile(r2Out);
        reads[2]->printFile(uOut);
        node->leaf->cluster->visited = true;
      }
    }
  }
  r1Out.close();
  r2Out.close();
  uOut.close();
  endMessage(log, start);

  for (Cluster* cluster: clusters) {
    delete cluster;
  }

  log
    << "\nRead " << line << " out of " << total << " lines of length "
      << length << " (" << (float)(total - line) / total << "% discarded).\n"
    << "Left after removing perfect duplicates: " << nonDuplicates << " ("
      << 100 * (float)nonDuplicates / line << "%).\n"
    << "Left after removing nonperfect duplicates (distance " << distance
      << "): " << id << " (" << 100 * (float)(id) / line << "%).\n";

  log.close();

}


/*
 * Argument parsing.
 */
int main(int argc, char* argv[]) {
  CliIO io(argc, argv);

  interface(
    io,
    dedup, argv[0], "Deduplicate a dataset.", 
      param("read1", "FastQ file for read 1"),
      param("read2", "FastQ file for read 2"),
      param("umi", "FastQ file for the UMI"),
      param("length", "word length"),
      param("r1", "FastQ file for read 1"),
      param("r2", "FastQ file for read 2"),
      param("u", "FastQ file for the UMI"),
      param("-d", 1, "distance"),
      param("-l", "/dev/stderr", "log file name"));

  return 0;
}
