#include <catch.hpp>

#include "../src/fastq.h"
#include "test_fastq.h"


TEST_CASE("Extract UMI from FastQ header") {
  // Read without an UMI in the header
  fakeRead *noUMI = new fakeRead(
    new string("@header"),
    new string("ATCG")
  );
  // Read with an UMI in the header
  fakeRead *hasUMI = new fakeRead(
    new string("@header_AATT"),
    new string("ATCG")
  );
  // Read with extra stuff after the UMI
  fakeRead *weirdHeader = new fakeRead(
    new string("@header_AATT with extra junk"),
    new string("ATCG")
  );
  // Read with multiple underscores
  fakeRead *underscoreHeader = new fakeRead(
    new string("@header_with_many_underscores_AATT and space"),
    new string("ATCG")
  );

  // Tests for read without UMI in header
  REQUIRE(_extractUMI(*noUMI->mName) == "");
  REQUIRE(!_hasUMI(*noUMI->mName));
  REQUIRE(_hasUMI(*hasUMI->mName));

  // Tests for reads with UMI in header
  REQUIRE(_extractUMI(*hasUMI->mName) == "AATT");
  REQUIRE(_extractUMI(*weirdHeader->mName) == "AATT");

  // Test for read with multiple underscores
  REQUIRE(_extractUMI(*underscoreHeader->mName) == "AATT");

  // Clean up
  delete noUMI;
  delete hasUMI;
  delete weirdHeader;
  delete underscoreHeader;
}
