#include <catch.hpp>

#include "../src/fastq.h"

string extractUMI_(string);
string makeStringSize_(string, size_t, char);


TEST_CASE("Extract UMI from header") {
  SECTION("Extract UMI from read pointer") {
    Read hasUMI("header_AATT", "", "", "", "");
    REQUIRE(extractUMI(&hasUMI) == "AATT");

    Read BCL("Instrument:RunID:FlowCellID:Lane:Tile:X:Y:ATCG", "", "", "");
    REQUIRE(extractUMI(&BCL) == "ATCG");
  }
  SECTION("Extract UMI with underscore") {
    SECTION("No UMI in header") {
      REQUIRE(extractUMI_("header") == "");
      REQUIRE(extractUMI_("header with spaces") == "");
      REQUIRE(extractUMI_("header_with_many_underscores and space") == "");
      REQUIRE(extractUMI_("header_ignore_lowercase_umi_aatt") == "");
      REQUIRE(extractUMI_("header space then_underscore") == "");
      REQUIRE(extractUMI_("header space then_underscore_AATT") == "");
    }

    SECTION("UMI in header") {
      REQUIRE(extractUMI_("header_AATT") == "AATT");
      REQUIRE(extractUMI_("header_AATT with spaces") == "AATT");
      REQUIRE(extractUMI_("header_with_many_underscores_AATT") == "AATT");
      REQUIRE(extractUMI_("header_with_many_underscores_AATT and space") == "AATT");
    }
  }

  SECTION("Extract UMI with colon") {
    SECTION("No UMI in header") {
      REQUIRE(extractUMI_("Instrument:RunID:FlowCellID:Lane:Tile:X:Y more stuf") == "");
    }

    SECTION("UMI in header") {
      REQUIRE(extractUMI_("Instrument:RunID:FlowCellID:Lane:Tile:X:Y:ATCG") == "ATCG");
      REQUIRE(extractUMI_("Instrument:RunID:FlowCellID:Lane:Tile:X:Y:ATCG more stuf") == "ATCG");
      REQUIRE(extractUMI_("Instrument:RunID:FlowCellID:Lane:Tile:X:Y:ATCG more_underscore") == "ATCG");
    }
  }
}

TEST_CASE("Test making a Word out of a vector of Reads") {
  Read read1("header", "AAAA", "", "");
  Read read2("header2", "TTTT", "", "");
  vector<Read*> reads { &read1, &read2 };

  Word word = makeWord(reads, {4, 4}, 0);
  vector<uint8_t> expected = { 0, 0, 0, 0, 3, 3, 3, 3};
  REQUIRE(word.data == expected);
}

TEST_CASE("Test padding) when fetching more than header UMI length") {
  // Test reads
  Read read1("header_AAAA", "TTTT", "", "");
  Read read2("header", "GGGG", "", "");
  vector<Read*> reads{ &read1, &read2 };

  // To store the results
  vector<char> nuc;
  string nucleotides;
  string expected;

  SECTION("Extract full reads and UMI from the header") {
    nuc = getNucleotides(reads, {4, 4}, 4);
    nucleotides = string(nuc.data(), nuc.size());
    expected = "AAAATTTTGGGG";
    REQUIRE(nucleotides == expected);
  }

  SECTION("Test padding of UMI header") {
    nuc = getNucleotides(reads, {4, 4}, 6);
    nucleotides = string(nuc.data(), nuc.size());
    expected = "AAAANNTTTTGGGG";
    REQUIRE(nucleotides == expected);
  }

  SECTION("Test padding of regular reads") {
    nuc = getNucleotides(reads, {5, 5}, 4);
    nucleotides = string(nuc.data(), nuc.size());
    expected = "AAAATTTTNGGGGN";
    REQUIRE(nucleotides == expected);
  }

  SECTION("Test extracting a subset of UMI from the header") {
    nuc = getNucleotides(reads, {0, 0}, 3);
    nucleotides = string(nuc.data(), nuc.size());
    expected = "AAA";
    REQUIRE(nucleotides == expected);
  }

  SECTION("Test extracting a subset from the reads") {
    nuc = getNucleotides(reads, {2, 2}, 0);
    nucleotides = string(nuc.data(), nuc.size());
    expected="TTGG";
    REQUIRE(nucleotides == expected);
  }

  SECTION("Test extracting unequal number of nucleotides from the reads") {
    nuc = getNucleotides(reads, {1, 3}, 0);
    nucleotides = string(nuc.data(), nuc.size());
    expected="TGGG";
    REQUIRE(nucleotides == expected);
  }
}

TEST_CASE("Test dividing nucleotides over files") {
  // Vector for the expected result
  vector<size_t> expected;

  //1 file, 10 nt
  expected = { 10 };
  REQUIRE(ntFromFile(1, 10) == expected);

  //3 files, 1 nt
  expected = { 0, 0, 1 };
  REQUIRE(ntFromFile(3, 1) == expected);

  //3 files, 2 nt
  expected = { 0, 0, 2 };
  REQUIRE(ntFromFile(3, 2) == expected);

  //3 files, 3 nt
  expected = { 1, 1, 1 };
  REQUIRE(ntFromFile(3, 3) == expected);

  //3 files, 13 nt
  expected = { 4, 4, 5 };
  REQUIRE(ntFromFile(3, 13) == expected);

  //3 files, 12 nt
  expected = { 4, 4, 4 };
  REQUIRE(ntFromFile(3, 12) == expected);

  //3 files, 11 nt
  expected = { 3, 3, 5 };
  REQUIRE(ntFromFile(3, 11) == expected);

  //3 files, 10 nt
  expected = { 3, 3, 4 };
  REQUIRE(ntFromFile(3, 10) == expected);

  //3 files, 9 nt
  expected = { 3, 3, 3 };
  REQUIRE(ntFromFile(3, 9) == expected);
}

TEST_CASE("Test extracting only the large UMI from the header"){
  Read read("header_AAAAAA", "TTTT", "", "");
  vector<Read*> reads{&read};

  vector<char> nuc = getNucleotides(reads, {0}, 4);
  string nucleotides = string(nuc.data(), nuc.size());
  string expected = "AAAA";

  REQUIRE(nucleotides == expected);
}

TEST_CASE("Test if a string is a valid UMI") {
  SECTION("Invalid UMIs") {
    REQUIRE(not validUMI(""));
    REQUIRE(not validUMI("atcg"));
    REQUIRE(not validUMI("ATCGP"));
    REQUIRE(not validUMI("1234"));
    REQUIRE(not validUMI("ATCGN"));
  }

  SECTION("Valid UMIs") {
    REQUIRE(validUMI("A"));
    REQUIRE(validUMI("ATCG"));
  }
}


TEST_CASE("Test extracting the last field from a string") {
  SECTION("Missing last field") {
    REQUIRE(extractLastField("", ':') == "");
    REQUIRE(extractLastField("nothing", ':') == "");
    REQUIRE(extractLastField("empty:", ':') == "");
  }

  SECTION("Including last field") {
    REQUIRE(extractLastField("last:field", ':') == "field");
    REQUIRE(extractLastField("three:differient:fields", ':') == "fields");
  }
}

TEST_CASE("Test making a string a given size") {
  REQUIRE(makeStringSize_("AA", 0, 'N') == "");
  REQUIRE(makeStringSize_("AA", 1, 'N') == "A");
  REQUIRE(makeStringSize_("AA", 2, 'N') == "AA");
  REQUIRE(makeStringSize_("AA", 3, 'N') == "AAN");
}
