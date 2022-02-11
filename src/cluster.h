#pragma once

#include <map>
#include <vector>

#include <stdlib.h>

using std::map;
using std::vector;

struct NLeaf;

/*!
 * Cluster structure.
 */
struct Cluster {
  size_t id;
  bool visited = false;

  size_t maxCount = 0;
  NLeaf* maxLeaf = NULL;
  size_t size = 0;

  Cluster(size_t);
};


void assignCluster(NLeaf*, Cluster*);
map<size_t, size_t> clusterStats(vector<Cluster*>&);
void freeClusters(vector<Cluster*>&);
