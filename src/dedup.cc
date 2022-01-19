#include <string>

#include "../lib/commandIO/src/commandIO.h"

#include "../../src/trie.tcc"
#include "../lib/ngs/ngs.h"

using std::ios;
using std::ofstream;
using std::string;

struct NLeaf : Leaf {
  vector<size_t> neighbours;
};


void logMessage(ofstream& log, char const message[]) {
  log << time(NULL) << ' ' << message;
  log.flush();
}

void dedup(
    string inputName, size_t length, string outputName, string logName) {
  Trie<4, NLeaf> trie;

  ofstream log(logName.c_str(), ios::in | ios::binary);
  logMessage(log, "Reading data.\n");
  size_t line = 0;
  for (vector<uint8_t> word: readFile(inputName.c_str(), length)) {
    NLeaf* leaf = trie.add(word);
    leaf->neighbours.push_back(line++);
  }

  logMessage(log, "Calculating neighbours.\n");
  for (Result<NLeaf> walkResult: trie.walk()) {
    for (Result<NLeaf> hammingResult: trie.hamming(walkResult.path, 1)) {
      if (walkResult.leaf != hammingResult.leaf) {
        walkResult.leaf->neighbours.insert(
          walkResult.leaf->neighbours.end(),
          hammingResult.leaf->neighbours.begin(),
          hammingResult.leaf->neighbours.end());
      }
    }
  }

  logMessage(log, "Writing results.\n");
  ofstream output(outputName.c_str(), ios::in | ios::binary);
  for (Result<NLeaf> result: trie.walk()) {
    output << result.leaf->neighbours[0] << ':';
    for (size_t neighbour: result.leaf->neighbours) {
      output << ' ' << neighbour;
    }
    output << '\n';
  }
  output.close();

  logMessage(log, "Done.\n");
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
