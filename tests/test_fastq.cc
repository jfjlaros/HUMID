#include <catch.hpp>

#include "../src/fastq.h"
#include "test_fastq.h"


TEST_CASE("Dummy test") {
  fakeRead *read = new fakeRead(
    new string("@header"),
    new string("ATCG")
  );

  // Access 'read' properties through the pointer, to mimic usage of fastp
  REQUIRE(*read->mName == "@header");
  REQUIRE(*read->mSeq == "ATCG");
  delete read;
}

TEST_CASE("Test if FastQ header contains a UMI") {
  fakeRead *noUMI = new fakeRead(
    new string("@header"),
    new string("ATCG")
  );

  fakeRead *hasUMI = new fakeRead(
    new string("@header_ATTT"),
    new string("ATCG")
  );


  REQUIRE(!_hasUMI(*noUMI->mName));
  REQUIRE(_hasUMI(*hasUMI->mName));

  delete noUMI; delete hasUMI;
}
