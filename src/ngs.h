#ifndef NGS_H_
#define NGS_H_

#include <fstream>
#include <vector>

using std::ifstream;
using std::vector;


vector<uint8_t> readWord(ifstream&, size_t);

#endif
