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
  size_t id;
  size_t maxCount = 0;
  struct NLeaf* maxLeaf = nullptr;
  size_t size = 0;
  bool visited = false;

  Cluster(size_t);
};

void assignMaxCluster(NLeaf*, Cluster*);
void assignDirectionalCluster(NLeaf*, Cluster*);
map<size_t, size_t> clusterStats(vector<Cluster*>&);
void freeClusters(vector<Cluster*>&);
