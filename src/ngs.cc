#include <iostream>
#include <map>

#include "ngs.h"

using std::cout;
using std::ios;
using std::map;

map<char, uint8_t> nuc = {{'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}};


vector<Read> read(vector<FastqReader*>& readers) {
  vector<Read> reads;
  for (FastqReader* reader: readers) {
    reads.push_back(*reader->read());
  }
  return reads;
}

/*!
 * Read a word from multiple files.
 *
 * \param readers FastQ readers.
 * \param size Word length.
 *
 * \return Word.
 */
Word read(vector<FastqReader*>& readers, size_t length) {
  Word result;
  for (FastqReader* reader: readers) {
    Read* read = reader->read();
    if (!read) {
      result.data.clear();
      return result;
    }
    for (size_t i = 0; i < length; i++) {
      char nucleotide = (*read->mSeq)[i];
      if (nuc.contains(nucleotide)) {
        result.data.push_back(nuc[nucleotide]);
      }
      else {
        result.data.push_back(nuc['G']);
        result.filtered = true;
      }
    }
    delete read;
  }
  return result;
}

/*!
 * Read all words from multiple files.
 *
 * \param read1 File name for read 1.
 * \param read2 File name for read 2.
 * \param umi File name for the UMI.
 * \param size Word length.
 *
 * \return All words.
 */
generator<Word> readFiles(
    string read1, string read2, string umi, size_t length) {
  FastqReader read1Reader(read1.c_str());
  FastqReader read2Reader(read2.c_str());
  FastqReader umiReader(umi.c_str());
  vector<FastqReader*> readers = {&read1Reader, &read2Reader, &umiReader};

  Word result = read(readers, length);
  while (result.data.size()) {
    co_yield result;
    result = read(readers, length);
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
