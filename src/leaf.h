#pragma once

#include "../lib/trie/src/trie.tcc"

/*!
 * Leaf structure for neighbour finding.
 */
struct NLeaf : Leaf {
  vector<size_t> lines;
  vector<NLeaf*> neighbours;
  Cluster* cluster = NULL;
};
