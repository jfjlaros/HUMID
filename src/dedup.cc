#include <string>

#include "../lib/commandIO/src/commandIO.h"

#include "../../src/trie.tcc"
#include "../lib/ngs/ngs.h"

using std::ios;
using std::ofstream;
using std::string;

struct NLeaf : Leaf {
  vector<size_t> lines;
  vector<NLeaf*> neighbours;
  size_t cluster = 0;
};


void assignCluster(NLeaf* leaf, size_t cluster) {
  leaf->cluster = cluster;
  for (NLeaf* neighbour: leaf->neighbours) {
    if (!neighbour->cluster) {
      assignCluster(neighbour, cluster);
    }
  }
}

time_t startMessage(ofstream& log, char const message[]) {
  log << message << "... ";
  log.flush();

  return time(NULL);
}

void endMessage(ofstream& log, time_t start) {
  unsigned int seconds = (unsigned int)difftime(time(NULL), start);
  log << "done. (" << seconds / 60 << 'm' << seconds % 60 << "s)\n";
  log.flush();
}

void dedup(
    string inputName, size_t length, string outputName, string logName) {
  Trie<4, NLeaf> trie;

  ofstream log(logName.c_str(), ios::out | ios::binary);
  time_t start = startMessage(log, "Reading data");
  size_t line = 0;
  for (vector<uint8_t> word: readFile(inputName.c_str(), length)) {
    NLeaf* leaf = trie.add(word);
    leaf->lines.push_back(line++);
  }
  endMessage(log, start);

  start = startMessage(log, "Calculating neighbours");
  for (Result<NLeaf> walkResult: trie.walk()) {
    for (Result<NLeaf> hammingResult: trie.hamming(walkResult.path, 1)) {
      if (walkResult.leaf != hammingResult.leaf) {
        walkResult.leaf->neighbours.push_back(hammingResult.leaf);
        hammingResult.leaf->neighbours.push_back(walkResult.leaf);
      }
    }
  }
  endMessage(log, start);

  start = startMessage(log, "Calculating clusters");
  size_t cluster = 1;
  for (Result<NLeaf> result: trie.walk()) {
    if (!result.leaf->cluster) {
      assignCluster(result.leaf, cluster++);
    }
  }
  endMessage(log, start);

  startMessage(log, "Writing results");
  ofstream output(outputName.c_str(), ios::out | ios::binary);
  for (Result<NLeaf> result: trie.walk()) {
    for (size_t line: result.leaf->lines) {
      output << line << ' ' << result.leaf->cluster << '\n';
    }
  }
  output.close();
  endMessage(log, start);

  log.close();
}


int main(int argc, char* argv[]) {
  CliIO io(argc, argv);

  interface(
    io,
    dedup, argv[0], "Deduplicate a dataset.", 
      param("input", "input file name"),
      param("length", "word length"),
      param("-o", "/dev/stdout", "output file name"),
      param("-l", "/dev/stderr", "log file name"));

  return 0;
}
