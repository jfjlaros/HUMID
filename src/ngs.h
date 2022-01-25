#pragma once

#include <fstream>
#include <vector>

#include "../../../src/CPP20Coroutines/include/generator.hpp"
#include "../fastp/src/fastqreader.h"

using std::ifstream;
using std::vector;

/*!
 * Word.
 */
struct Word {
  vector<uint8_t> data;
  bool filtered = false;
};

generator<Word> readFiles(string, string, string, size_t);
void printWord(vector<uint8_t>&);
