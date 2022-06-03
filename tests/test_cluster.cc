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
