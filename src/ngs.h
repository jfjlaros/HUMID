#pragma once

#include <vector>

#include "../../../src/CPP20Coroutines/include/generator.hpp"
#include "../fastp/src/fastqreader.h"

using std::vector;

/*!
 * Word.
 */
struct Word {
  vector<uint8_t> data;
  bool filtered = false;
};

generator<vector<Read*>> readFiles(vector<string>&);
Word makeWord(vector<Read*>&, size_t);
void printWord(vector<uint8_t>&);
