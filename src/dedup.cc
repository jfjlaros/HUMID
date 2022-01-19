#include <iostream>

#include "../../src/trie.tcc"
#include "../lib/ngs/ngs.h"

using std::cerr;
using std::cout;

char const name[] = "../../data/100k.bin";

struct NLeaf : Leaf {
  vector<size_t> neighbours;
};


int main(void) {
  Trie<4, NLeaf> trie;

  cerr << time(NULL) << " Reading data.\n";
  size_t line = 0;
  for (vector<uint8_t> word: readFile(name, 24)) {
    NLeaf* leaf = trie.add(word);
    leaf->neighbours.push_back(line++);
  }

  cerr << time(NULL) << " Calculating neighbours.\n";
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

  cerr << time(NULL) << " Writing results.\n";
  for (Result<NLeaf> result: trie.walk()) {
    cout << result.leaf->neighbours[0] << ':';
    for (size_t neighbour: result.leaf->neighbours) {
      cout << ' ' << neighbour;
    }
    cout << '\n';
  }

  cerr << time(NULL) << " Done.\n";

  return 0;
}
