#include <catch.hpp>

#include "../src/cluster.h"
#include "../src/leaf.h"

bool atLeastDouble_(int, int);
bool atMostHalf_(int, int);
NLeaf* maxNeighbour(NLeaf*);


// Helper function to link nodes
void link(NLeaf* a, NLeaf* b) {
  a->neighbours.push_back(b);
  b->neighbours.push_back(a);
}


TEST_CASE("Test if a is at least 2x b", "[cluster]") {
  REQUIRE(atLeastDouble_(1, 0));
  REQUIRE(atLeastDouble_(2, 1));
  REQUIRE(not atLeastDouble_(3, 2));
}

TEST_CASE("Test if a is at most half of b", "[cluster]") {
  REQUIRE(atMostHalf_(0, 1));
  REQUIRE(atMostHalf_(1, 2));
  REQUIRE(not atMostHalf_(2, 3));
}

TEST_CASE("Test walking a node with no neighbours", "[cluster]") {
  // Create a node that is all alone
  NLeaf alone;
  //A leaf with no neighbours should return itself
  REQUIRE(maxNeighbour(&alone) == &alone);
}

TEST_CASE("Test walking node whose neighbour is already assigned", "[cluster]") {
  // Create a more complex setup, of chained leafs
  NLeaf leaf;
  leaf.count = 1;

  //A neighbour that is already in a cluster should not be used
  NLeaf assigned_neighbour;
  Cluster c(2);
  assigned_neighbour.cluster = &c;
  assigned_neighbour.count = 2;

  link(&leaf, &assigned_neighbour);
  REQUIRE(maxNeighbour(&leaf) == &leaf);
}

TEST_CASE("Test walking a chain of nodes", "[cluster]") {
  //If there is a neighbour that is not assigned and it conforms to the 2x
  //requirement, it should be used
  NLeaf leaf;
  leaf.count = 1;

  NLeaf free_neighbour;
  free_neighbour.count = 2;

  link(&leaf, &free_neighbour);
  REQUIRE(maxNeighbour(&leaf) == &free_neighbour);

  //Lets test a third, further neighbour
  NLeaf third_neighbour;
  third_neighbour.count = 4;
  link(&free_neighbour, &third_neighbour);
  REQUIRE(maxNeighbour(&leaf) == &third_neighbour);

  //Add one more neighbour, that is not high enough to add
  NLeaf last_one;
  last_one.count = 7;
  link(&third_neighbour, &last_one);

  // Check that the last neighbour was not added, since it was not high
  // enough
  REQUIRE(maxNeighbour(&leaf) == &third_neighbour);
}

TEST_CASE("Test assigning to cluster", "[cluster]") {
  // Initialise the node counts
  NLeaf node1; // smallest node
  node1.count = 2;

  NLeaf node2; // intermediate node
  node2.count = 4;

  NLeaf node3; // largest reachable node by walking up in 2x steps
  node3.count = 8;

  NLeaf node4; // largest overall node, only reachable from node5
  node4.count = 10;

  NLeaf node5; // Only reachable from node5
  node5.count = 3;

  // Test the inialisation
  REQUIRE(node4.count == 10);

  // Link node1 and node2 together
  link(&node1, &node2);

  // Test that node1 and node2 are now linked
  REQUIRE(node1.neighbours.front() == &node2);
  REQUIRE(node2.neighbours.front() == &node1);

  // Link the rest of the nodes
  link(&node2, &node3);
  link(&node3, &node4);
  link(&node4, &node5);

  // Assume we start with node 1, and want to assign them to cluster 1
  Cluster cluster1(1);

  // Next, we assign all nodes
  assignDirectionalCluster(&node1, &cluster1);

  // Then, the first three nodes should be assigned to cluster 1
  REQUIRE(node1.cluster->id == 1);
  REQUIRE(node2.cluster->id == 1);
  REQUIRE(node3.cluster->id == 1);

  // The last two node should not be assigned
  REQUIRE(not node4.cluster);
  REQUIRE(not node5.cluster);

  // Now we assign node4 to cluster2
  Cluster cluster2(2);
  assignDirectionalCluster(&node4, &cluster2);
  REQUIRE(node4.cluster->id == 2);
  REQUIRE(node5.cluster->id == 2);

  // Check that the cluster size is correct
  REQUIRE(cluster1.size == 14);
  REQUIRE(cluster2.size == 13);

  // Check that the maxleaf is correct
  REQUIRE(cluster1.maxLeaf == &node3);
  REQUIRE(cluster2.maxLeaf == &node4);

  // Check that the maxCount is correct
  REQUIRE(cluster1.maxCount == 8);
  REQUIRE(cluster2.maxCount == 10);
}
