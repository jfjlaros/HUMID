#pragma once

#include <map>
#include <vector>

#include <stdlib.h>

using std::map;
using std::vector;

/*!
 * Cluster structure.
 */
struct Cluster {
  size_t const id;
  size_t maxCount = 0;
  struct NLeaf* maxLeaf = nullptr;
  size_t size = 0;
  bool visited = false;

  Cluster(size_t const);
};

void assignMaxCluster(NLeaf*, Cluster* const);
void assignDirectionalCluster(NLeaf* const, Cluster*);
map<size_t, size_t> clusterStats(vector<Cluster*> const&);
void freeClusters(vector<Cluster*> const&);
