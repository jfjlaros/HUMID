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

  Word word = makeWord(reads, 8);
  vector<uint8_t> expected = { 0, 0, 0, 0, 3, 3, 3, 3};
  REQUIRE(word.data == expected);
}

TEST_CASE("Test extracting nucleotides from a vector of Reads") {
  Read read1("header", "AAAA", "", "");
  Read read2("header2", "TTTT", "", "");
  vector<Read*> reads { &read1, &read2};

  vector<char> nuc = getNucleotides(reads, 8);
  string nucleotides = string(nuc.data(), nuc.size());
  string expected = "AAAATTTT";

  REQUIRE(nucleotides == expected);
}

TEST_CASE("Test extracting UMI from read when UMI is longer than wordSize") {
  Read read("header_AAAA", "TTTT", "", "");
  vector<Read*> reads{ &read };

  vector<char> nuc = getNucleotides(reads, 3);
  string nucleotides = string(nuc.data(), nuc.size());
  string expected = "AAA";

  REQUIRE(nucleotides == expected);
}

TEST_CASE("Test extracting nucleotides from header and read"){
  Read read("header_AAAA", "TTTT", "", "");
  vector<Read*> reads{ &read };

  vector<char> nuc = getNucleotides(reads, 6);
  string nucleotides = string(nuc.data(), nuc.size());
  string expected = "AAAATT";

  REQUIRE(nucleotides == expected);
}

TEST_CASE("Test dividing nucleotides over files") {
  // Vector for the expected result
  vector<size_t> expected;

  //1 file, 10 nt
  expected = { 10 };
  REQUIRE(_ntFromFile(1, 10) == expected);

  //3 files, 1 nt
  expected = { 0, 0, 1 };
  REQUIRE(_ntFromFile(3, 1) == expected);

  //3 files, 2 nt
  expected = { 0, 0, 2 };
  REQUIRE(_ntFromFile(3, 2) == expected);

  //3 files, 3 nt
  expected = { 1, 1, 1 };
  REQUIRE(_ntFromFile(3, 3) == expected);

  //3 files, 13 nt
  expected = { 4, 4, 5 };
  REQUIRE(_ntFromFile(3, 13) == expected);

  //3 files, 12 nt
  expected = { 4, 4, 4 };
  REQUIRE(_ntFromFile(3, 12) == expected);

  //3 files, 11 nt
  expected = { 3, 3, 5 };
  REQUIRE(_ntFromFile(3, 11) == expected);

  //3 files, 10 nt
  expected = { 3, 3, 4 };
  REQUIRE(_ntFromFile(3, 10) == expected);

  //3 files, 9 nt
  expected = { 3, 3, 3 };
  REQUIRE(_ntFromFile(3, 9) == expected);
}

TEST_CASE("Test extracting nucleotides from header and read when wordSize has remainder"){
  Read read1("header_AAAA", "TTTT", "", "");
  Read read2("header", "GGGG", "", "");
  vector<Read*> reads{ &read1, &read2 };

  vector<char> nuc = getNucleotides(reads, 11);
  string nucleotides = string(nuc.data(), nuc.size());
  string expected = "AAAATTTGGGG";

  REQUIRE(nucleotides == expected);
}

TEST_CASE("Test extracting only the large UMI from the header"){
  Read read("header_AAAAAA", "TTTT", "", "");
  vector<Read*> reads{&read};

  vector<char> nuc = getNucleotides(reads, 4);
  string nucleotides = string(nuc.data(), nuc.size());
  string expected = "AAAA";

  REQUIRE(nucleotides == expected);
}

TEST_CASE("Test finding all occurrences of char in string") {
  vector<size_t> expected;

  expected = { 2, 5, 6, 9 };
  REQUIRE(findAll(':',"01:34::78:") == expected);
}

TEST_CASE("Test if a string is a valid UMI") {
  // Invalid UMIs
  REQUIRE(not _validUMI(""));
  REQUIRE(not _validUMI("atcg"));
  REQUIRE(not _validUMI("ATCGP"));
  REQUIRE(not _validUMI("1234"));

  // Valid UMIs
  REQUIRE(_validUMI("A"));
  REQUIRE(_validUMI("ATCGN"));
}
