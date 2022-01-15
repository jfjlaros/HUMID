#include "ngs.h"


vector<uint8_t> readWord(ifstream& handle, size_t size) {
  uint8_t buffer[size];
  handle.read((char*)buffer, size);
  return vector<uint8_t>(buffer, buffer + size);
}
