#include <catch.hpp>

#include "../src/fastq.h"


TEST_CASE("Test extracting a UMI from a fastp Read pointer") {
  Read hasUMI("header_AATT", "", "", "", "");
  REQUIRE(extractUMI(&hasUMI) == "AATT");
}

TEST_CASE("Test extracting UMI from FastQ header") {
  // Tests for read without UMI in header
  REQUIRE(_extractUMI("header") == "");
  REQUIRE(_extractUMI("header with spaces") == "");
  REQUIRE(_extractUMI("header_with_many_underscores and space") == "");
  REQUIRE(_extractUMI("header_ignore_lowercase_umi_aatt") == "");
  REQUIRE(_extractUMI("header space then_underscore") == "");
  REQUIRE(_extractUMI("header space then_underscore_AATT") == "");

  // Tests for reads with UMI in header
  REQUIRE(_extractUMI("header_AATT") == "AATT");
  REQUIRE(_extractUMI("header_AATT with spaces") == "AATT");
  REQUIRE(_extractUMI("header_with_many_underscores_AATT") == "AATT");
  REQUIRE(_extractUMI("header_with_many_underscores_AATT and space") == "AATT");
}

TEST_CASE("Test making a Word out of a vector of Reads") {
  Read read1("header", "AAAA", "", "");
  Read read2("header2", "TTTT", "", "");
  vector<Read*> reads;
  reads.push_back(&read1);
  reads.push_back(&read2);
  REQUIRE(true);

  Word word = makeWord(reads, 4);
  vector<uint8_t> expected = { 0, 0, 0, 0, 3, 3, 3, 3};
  REQUIRE(word.data == expected);
}
