#pragma once

#include <string>
#include <vector>

#include "../lib/trie/lib/CPP20Coroutines/include/generator.hpp"
#include "../lib/fastp/src/fastqreader.h"

using std::string;
using std::vector;

/*!
 * Word.
 */
struct Word {
  vector<uint8_t> data {};
  bool filtered {false};
};

generator<vector<Read*>> readFiles(vector<string> const);
vector<char> getNucleotides(
  vector<Read*> const&, vector<size_t> const, size_t const);
Word makeWord(vector<Read*> const&, vector<size_t> const, size_t const);
void printWord(vector<uint8_t> const&);
string addDir(char const[], string const);
string makeFileName(string const, string const, string const);
vector<string> makeFileNames(vector<string> const, string const, string const);
string extractLastField(string const, char const);
bool validUMI(string const);
string extractUMI(Read* const);
vector<size_t> ntFromFile(size_t const, size_t const);
