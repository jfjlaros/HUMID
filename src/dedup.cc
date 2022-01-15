#include <fstream>
#include <ios>
#include <iostream>
#include <vector>

#include "../../src/trie.tcc"
#include "../lib/ngs/ngs.h"

using std::cout;
using std::ios;

const char file[] = "../../data/100k.bin";


class NLeaf : public Leaf<size_t> {
  public:
    size_t line = 0;

    void add(size_t& l) { line = l; }
};

void visit(vector<uint8_t>& word, NLeaf& leaf) {
  for (uint8_t letter: word) {
    cout << (int)letter << ' ';
  }
  cout << leaf.line << '\n';
}


int main(void) {
  Trie<4, NLeaf> trie;

  ifstream handle(file, ios::in | ios::binary);
  size_t line = 0;
  while (!handle.eof()) {
    vector<uint8_t> word = readWord(handle, 24);
    trie.add(word, line);
    line++;
    handle.peek();
  }
  handle.close();

  trie.traverse(visit);

  return 0;
}
