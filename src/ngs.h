#pragma once

#include <fstream>
#include <vector>

#include "../../../src/CPP20Coroutines/include/generator.hpp"
#include "../fastp/src/fastqreader.h"

using std::ifstream;
using std::vector;


vector<uint8_t> readWord(ifstream&, size_t);
vector<uint8_t> readWord(FastqReader[3], size_t);
generator<vector<uint8_t>> readFile(char const[], size_t);
generator<vector<uint8_t>> readFiles(string, string, string, size_t);
void printWord(vector<uint8_t>&);
