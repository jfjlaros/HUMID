#pragma once

#include <fstream>
#include <vector>

#include "../../../src/CPP20Coroutines/include/generator.hpp"

using std::ifstream;
using std::vector;


vector<uint8_t> readWord(ifstream&, size_t);
generator<vector<uint8_t>> readFile(char const[], size_t);
void printWord(vector<uint8_t>&);
