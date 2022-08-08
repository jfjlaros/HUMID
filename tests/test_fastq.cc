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
