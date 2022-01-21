#include <iostream>

#include "ngs.h"

using std::cout;
using std::ios;


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
