#include <iostream>
#include <map>
#include <unordered_set>

#include "fastq.h"

using std::cout;
using std::ios;
using std::map;
using std::unordered_set;

map<char, uint8_t> nuc = {{'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}};
unordered_set<char> nnuc = {'A', 'T', 'C', 'G', 'N'};

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
    if (not read) {
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
  while (not readVector.eof) {
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
 * Extract `wordLength` nucleotides from `reads`
 * If the first file has a UMI in the header, this will get preference.
 */

vector<char> getNucleotides(vector<Read*>& reads, size_t wordLength) {
  vector<char> nucleotides;

  // Pull the UMI from the header of the first read
  string headerUMI = extractUMI(reads.front());
  for (size_t i=0; i < wordLength and i < headerUMI.size(); i++) {
      nucleotides.push_back(headerUMI[i]);
    }

  // The length we still have available from wordLength after extracting the
  // UMI from the header
  size_t length;

  // If the UMI in the header is longer (or the same size) as the wordLength,
  // we are done
  if (wordLength <= headerUMI.size()) {
    return nucleotides;
  }
  else {
    length = (wordLength - headerUMI.size())/ reads.size();
  }

  //size_t length = wordLength / reads.size();
  for (Read* read: reads) {
    for (size_t i = 0; i < length; i++) {
      char nucleotide = (*read->mSeq)[i];
      nucleotides.push_back(nucleotide);
    }
  }
  return nucleotides;
}

/*!
 * Select a total of `wordLength` nucleotides from every read in `reads` to
 * create a word.
 *
 * \param reads Reads.
 * \param wordLength Read selection length.
 *
 * \return Word.
 */
Word makeWord(vector<Read*>& reads, size_t wordLength) {
  Word word;
  vector<char> nucleotides = getNucleotides(reads, wordLength);
  for (char nucleotide: nucleotides) {
    if (nuc.contains(nucleotide)) {
      word.data.push_back(nuc[nucleotide]);
    }
    else {
      word.data.push_back(nuc['G']);
      word.filtered = true;
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
 * Extract UMI from a header
 *
 * \param header Fastq header line
 */
string _extractUMI(string header) {
  size_t first_space = header.find(" ");
  size_t umiStart = header.substr(0, first_space).find_last_of("_");

  // We need to do some sanity checking on the eventual 'UMI' we find.
  string UMI;

  // If there is no underscore in the header, there is no UMI.
  if (umiStart == string::npos) {
    UMI="";
  }
  // Get the string between the last _ and the first space.
  else {
    UMI = header.substr(umiStart + 1, first_space - umiStart - 1);
  }
  // Check if the UMI only contains characters from 'ATCGN'. If we find any
  // other character, the 'UMI' is not actually a UMI.
  for (char c: UMI) {
    if (nnuc.find(c) == nnuc.end()) {
      return "";
    }
  }
  return UMI;
}

/*!
 * Extract UMI from a read
 *
 * \param read Read
 */
string extractUMI(Read* read) {
  return _extractUMI(*read->mName);
}
