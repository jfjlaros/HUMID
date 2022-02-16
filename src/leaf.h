#pragma once

#include "../lib/trie/src/trie.tcc"

/*!
 * Leaf structure for neighbour finding.
 */
struct NLeaf : Leaf {
  vector<NLeaf*> neighbours;
  Cluster* cluster = NULL;
};
