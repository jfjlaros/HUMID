#include <catch.hpp>
#include "../src/cluster.cc"

TEST_CASE("Test PCR step size") {
    NLeaf leaf;
    NLeaf neighbour1;
    NLeaf neighbour2;

    // Check that NLeaf count is initialised to 0
    REQUIRE(leaf.count == 0);

    // Set the counts
    leaf.count = 4;
    neighbour1.count = 2;

    REQUIRE(_pcr_step_size(&leaf, &neighbour1));
}

TEST_CASE("Test if a is at least 2x b", "[cluster]"){
    REQUIRE(_at_least_double(1, 0));
    REQUIRE(_at_least_double(2, 1));
    REQUIRE(!_at_least_double(3,2));
}

TEST_CASE("Test if a is at most half of b", "[cluster]"){
    REQUIRE(_at_most_half(0, 1));
    REQUIRE(_at_most_half(1, 2));
    REQUIRE(!_at_most_half(2, 3));
}

TEST_CASE("Test walking neighbours to find a maximum node", "[cluster]"){
    NLeaf leaf;
    leaf.count = 1;

    //A leaf with no neighbours should return itself
    //NLeaf* lp = &leaf;
    REQUIRE(_max_neighbour(&leaf) == &leaf);

    //A neighbour that is already in a cluster should not be used
    NLeaf assigned_neighbour;
    Cluster* c = new Cluster(2);
    assigned_neighbour.cluster = c;
    leaf.neighbours.push_back(&assigned_neighbour);
    REQUIRE(_max_neighbour(&leaf) == &leaf);

    //If there is a neighbour that is not assigned and it conforms to the 2x
    //requirement, it should be used
    NLeaf free_neighbour;
    free_neighbour.count=2;
    leaf.neighbours.push_back(&free_neighbour);
    REQUIRE(_max_neighbour(&leaf) == &free_neighbour);

    //Lets test a third, further neighbour
    NLeaf third_neighbour;
    third_neighbour.count=4;
    free_neighbour.neighbours.push_back(&third_neighbour);
    REQUIRE(_max_neighbour(&leaf) == &third_neighbour);

    //Add one more neighbour, that is not high enough to add
    NLeaf last_one;
    last_one.count=7;
    third_neighbour.neighbours.push_back(&last_one);
    // Should not be the last one, but the third one
    REQUIRE(_max_neighbour(&leaf) == &third_neighbour);
}
