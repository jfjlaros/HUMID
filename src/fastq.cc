#include <iostream>
#include <map>

#include "fastq.h"

using std::cout;
using std::ios;
using std::map;

map<char, uint8_t> nuc = {{'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}};

/*
 * Read vector.
 */
struct _ReadVector {
  vector<Read*> reads;
  bool eof = false;

  void free();
};


/*
 * Destroy a ReadVector.
 */
void _ReadVector::free() {
  for (Read* read: reads) {
    delete read;
  }
}


/*
 * Read one FastQ record from multiple files.
 *
 * \param readers FastQ readers.
 *
 * \return FastQ records plus EOF status.
 */
_ReadVector _readFastq(vector<FastqReader*>& readers) {
  _ReadVector readVector;
  for (FastqReader* reader: readers) {
    Read* read = reader->read();
    if (!read) {
      readVector.eof = true;
    }
    readVector.reads.push_back(read);
  }
  return readVector;
}

/*!
 * Loop over all reads in multiple FastQ files.
 *
 * \param files FastQ file names.
 *
 * \return All reads.
 */
generator<vector<Read*>> readFiles(vector<string> files) {
  vector<FastqReader*> readers;
  for (string file: files) {
    FastqReader* reader = new FastqReader(file.c_str());
    readers.push_back(reader);
  }

  _ReadVector readVector = _readFastq(readers);
  while (!readVector.eof) {
    co_yield readVector.reads;
    readVector.free();
    readVector = _readFastq(readers);
  }
  // TODO: Extra readVector.free() here when files are not of equal length?

  for (FastqReader* reader: readers) {
    delete reader;
  }
}

/*!
 * Select `length` nucleotides from every read in `reads` to create a word.
 *
 * \param reads Reads.
 * \param length Read selection length.
 *
 * \return Word.
 */
Word makeWord(vector<Read*>& reads, size_t length) {
  Word word;
  for (Read* read: reads) {
    for (size_t i = 0; i < length; i++) {
      char nucleotide = (*read->mSeq)[i];
      if (nuc.contains(nucleotide)) {
        word.data.push_back(nuc[nucleotide]);
      }
      else {
        word.data.push_back(nuc['G']);
        word.filtered = true;
      }
    }
  }
  return word;
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

/*!
 */
string addDir(char const filename[], string dir) {
  return dir + '/' + filename;
}

/*!
 */
string makeFileName(string filename, string dir, string suffix) {
  string name = basename(filename.c_str());
  size_t pos = name.find('.');
  string suff =
    name.substr(0, pos) + '_' + suffix + name.substr(pos, string::npos);
  return addDir(suff.c_str(), dir);
}

/*!
 */
vector<string> makeFileNames(vector<string> files, string dir, string suffix) {
  vector<string> fileNames;
  for (string name: files) {
    fileNames.push_back(makeFileName(name, dir, suffix));
  }
  return fileNames;
}

/*!
 * Determine whether header contains a UMI before the first space
 *
 * \param header FastQ header line
 */
bool _hasUMI(string header) {
  return !_extractUMI(header).empty();
}

/*!
 * Extract UMI from a header
 *
 * \param header Fastq header line
 */
string _extractUMI(string header) {
  size_t first_space = header.find(" ");
  size_t umiStart = header.substr(0, first_space).find("_");

  // If there is no underscore in the header
  if (umiStart == string::npos) {
    return "";
  }
  else {
    return header.substr(umiStart + 1, first_space - umiStart - 1);
  }
}
