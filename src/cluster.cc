#include "cluster.h"
#include "leaf.h"


/*!
 * Constructor.
 */
Cluster::Cluster(size_t id) {
  this->id = id;
}


/*!
 * Assign leaf to cluster, and update cluster->maxLeaf if needed
 *
 * \param leaf Leaf node.
 * \param cluster Cluster.
 */
void _assignLeaf(NLeaf* leaf, Cluster* cluster) {
  leaf->cluster = cluster;
  leaf->cluster->size += leaf->count;
  if (leaf->count > cluster->maxCount) {
    cluster->maxLeaf = leaf;
    cluster->maxCount = leaf->count;
  }
}

/*!
 * Traverse neighbours to assign cluster IDs
 *
 * \param leaf Leaf node.
 * \param cluster Cluster.
 */
void assignMaxCluster(NLeaf* leaf, Cluster* cluster) {
  _assignLeaf(leaf, cluster);
  for (NLeaf* neighbour: leaf->neighbours) {
    if (!neighbour->cluster) {
      assignMaxCluster(neighbour, cluster);
    }
  }
}

/*!
 * Determine if the count difference between leaf and neighbour is indicative
 * of an error that has been amplified by PCR
 *
 * \param leaf Leaf node.
 * \param neighbour Leaf node.
 */
bool _pcr_step_size(NLeaf* leaf, NLeaf* neighbour) {
    return leaf->count > 2 * neighbour->count - 1;
}

/*!
 * Traverse neighbours to assign cluster IDs, using the directional method
 * Only add neighbours that have at least 2x more reads
 *
 * \param leaf Leaf node.
 * \param cluster Cluster.
 */
void _assignDirectionalClusterUp(NLeaf* leaf, Cluster* cluster){
  _assignLeaf(leaf, cluster);

  for (NLeaf* neighbour: leaf->neighbours) {
    if (!neighbour->cluster) {
      // If the neighbour has more than 2x the number of reads, the neighbour
      // is the 'true' sequence
      if (_pcr_step_size(neighbour, leaf))
        _assignDirectionalClusterUp(neighbour, cluster);
    }
  }
}

/*!
 * Traverse neighbours to assign cluster IDs, using the directional method
 * Only add neighbours that have at most 2x fewer reads
 *
 * \param leaf Leaf node.
 * \param cluster Cluster.
 */
void _assignDirectionalClusterDown(NLeaf* leaf, Cluster* cluster){
  _assignLeaf(leaf, cluster);

  for (NLeaf* neighbour: leaf->neighbours) {
    if (!neighbour->cluster) {
      // If the neighbour has more than 2x the number of reads, the neighbour
      // is the 'true' sequence
      if (_pcr_step_size(leaf, neighbour))
        _assignDirectionalClusterDown(neighbour, cluster);
    }
  }
}

/*!
 * Traverse neighbours to assign cluster IDs, using the directional method
 *
 * \param leaf Leaf node.
 * \param cluster Cluster.
 */
void assignDirectionalCluster(NLeaf* leaf, Cluster* cluster){
  _assignLeaf(leaf, cluster);

  for (NLeaf* neighbour: leaf->neighbours) {
    if (!neighbour->cluster) {
      // If the neighbour has less than half of the number of reads, the
      // neighbour belongs to the current cluster
      if (_pcr_step_size(leaf, neighbour)) {
        _assignDirectionalClusterDown(neighbour, cluster);
      }
      // If the neighbour has more than 2x the number of reads, the neighbour
      // is the 'true' sequence
      else if (_pcr_step_size(neighbour, leaf)) {
        _assignDirectionalClusterUp(neighbour, cluster);
      }
    }
  }
}

/*!
 * Make a histogram of cluster sizes.
 *
 * \param clusters List of clusters.
 *
 * \return Histogram of cluster sizes.
 */
map<size_t, size_t> clusterStats(vector<Cluster*>& clusters) {
  map<size_t, size_t> counts;
  for (Cluster* cluster: clusters) {
    counts[cluster->size]++;
  }
  return counts;
}

/*!
 * Destroy a list of clusters.
 *
 * \param clusters List of clusters.
 */
void freeClusters(vector<Cluster*>& clusters) {
  for (Cluster* cluster: clusters) {
    delete cluster;
  }
}
