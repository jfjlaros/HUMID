#include <iostream>
#include <map>

#include "ngs.h"

using std::cout;
using std::ios;
using std::map;

map<char, uint8_t> nuc = {{'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}};

/*!
 * Read a word from a file.
 *
 * \param handle File handle.
 * \param size Word length.
 *
 * \return Word.
 */
vector<uint8_t> readWord(ifstream& handle, size_t size) {
  uint8_t buffer[size];
  handle.read((char*)buffer, size);
  return vector<uint8_t>(buffer, buffer + size);
}

/*!
 * Read a word from three files.
 *
 * \param result Word.
 * \param readers FastQ readers.
 * \param size Word length.
 *
 * \return -1: EOF, 0: success, 1: invalid read.
 */
int _filteredReadWord(
    vector<uint8_t>& result, FastqReader readers[3], size_t length) {
  for (uint8_t i = 0; i < 3; i++) {
    Read *read = readers[i].read();
    if (!read) {
      return -1;
    }
    for (size_t j = 0; j < length; j++) {
      if (nuc.contains((*read->mSeq)[j])) {
        result.push_back(nuc[(*read->mSeq)[j]]);
      }
      else {
        return 1;
      }
    }
  }
  return 0;
}

/*!
 * Read a word from three files.
 *
 * \param readers FastQ readers.
 * \param size Word length.
 *
 * \return Word.
 */
vector<uint8_t> readWord(FastqReader readers[3], size_t length) {
  vector<uint8_t> result;

  int status = 0;
  while (status > -1) {
    status = _filteredReadWord(result, readers, length);
    if (!status) {
      return result;
    }
    result.clear();
  }

  return vector<uint8_t>();
}

/*!
 * Read all words in a file.
 *
 * \param name File name.
 * \param size Word length.
 *
 * \return All words.
 */
generator<vector<uint8_t>> readFile(char const name[], size_t size) {
  ifstream handle(name, ios::in | ios::binary);
  while (!handle.eof()) {
    co_yield readWord(handle, 24);
    handle.peek();
  }
  handle.close();
}

/*!
 * Read all words in three files.
 *
 * \param read1 File name for read 1.
 * \param read2 File name for read 2.
 * \param umi File name for the UMI.
 * \param size Word length.
 *
 * \return All words.
 */
generator<vector<uint8_t>> readFiles(
    string read1, string read2, string umi, size_t length) {
  FastqReader readers[3] = {
    FastqReader(read1.c_str()), FastqReader(read2.c_str()),
    FastqReader(umi.c_str())};

  vector<uint8_t> result = readWord(readers, length);
  while (result.size()) {
    co_yield result;
    result = readWord(readers, length);
  }
}

/*!
 * Print a word.
 *
 * \param word Word.
 */
void printWord(vector<uint8_t>& word) {
  for (uint8_t letter: word) {
    cout << ' ' << (int)letter;
  }
  cout << '\n';
}
