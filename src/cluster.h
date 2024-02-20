#pragma once

#include <map>
#include <vector>

#include <stdlib.h>

using std::map;
using std::vector;

/*! Cluster structure. */
struct Cluster {
  size_t id;
  size_t maxCount {0};
  struct NLeaf* maxLeaf {nullptr};
  size_t size {0};
  bool visited {false};
};


/*! Traverse neighbours to assign cluster IDs.
 *
 * \param leaf Leaf node.
 * \param cluster Cluster.
 */
void assignMaxCluster(NLeaf*, Cluster* const);

/*! Traverse neighbours to assign cluster IDs, using the directional method.
 *
 * Also updates the maxCount for the cluster after determining a maximum
 * neighbour.
 *
 * \param leaf Leaf node.
 * \param cluster Cluster.
 */
void assignDirectionalCluster(NLeaf* const, Cluster*);

/*! Make a histogram of cluster sizes.
 *
 * \param clusters List of clusters.
 *
 * \return Histogram of cluster sizes.
 */
map<size_t, size_t> clusterStats(vector<Cluster*> const&);

/*! Destroy a list of clusters.
 *
 * \param clusters List of clusters.
 */
void freeClusters(vector<Cluster*> const&);
