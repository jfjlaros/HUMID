#include <iostream>
#include <map>
#include <sstream>

#include "fastq.h"
#include "../lib/fastp/src/util.h"

using std::cout;
using std::ios;
using std::map;

map<char const, uint8_t const> nuc {{'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}};

/*
 * Read vector.
 */
struct ReadVector_ {
  void free() const;

  vector<Read*> reads {};
  bool eof {false};
};


/*
 * Destroy a ReadVector.
 */
void ReadVector_::free() const {
  for (Read* const read: reads) {
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
ReadVector_ readFastq_(vector<FastqReader*> const& readers) {
  ReadVector_ readVector;
  for (FastqReader* const reader: readers) {
    Read* read {reader->read()};
    if (not read) {
      readVector.eof = true;
    }
    readVector.reads.push_back(read);
  }
  return readVector;
}

/*
 * Make string s the specified size, by either cutting it, or padding it.
 *
 * \param s String to make a certain size.
 * \param size Size the specified string should be.
 * \param padding Character to use for padding.
 *
 * \return New string made to size s.
 */
string makeStringSize_(string s, size_t const size, char const padding) {
  if (size == s.size()) {
    return s;
  }
  if (size < s.size()) {
    return s.substr(0, size);
  }
  // Add padding.
  return s.append(size - s.size(), padding);
}

/*
 * Extract UMI from a header.
 *
 * \param header Fastq header line.
 */
string extractUMI_(string const header) {
  size_t first_space {header.find(" ")};

  // The UMI must be before the first space.
  string substr {header.substr(0, first_space)};

  // If we detect a UMI with a _ separator, we return that UMI.
  string umi {extractLastField(substr, '_')};

  // Check if the UMI only contains characters from 'ATCGN'. If we find any
  // other character, the 'UMI' is not actually a UMI.
  if (validUMI(umi)) {
    return umi;
  }

  // Otherwise, we check to see if we find a BCL Convert style UMI.
  umi = extractLastField(substr, ':');
  if (validUMI(umi)) {
    return umi;
  }
  return "";
}


generator<vector<Read*>> readFiles(vector<string> const files) {
  vector<FastqReader*> readers;
  for (string const& file: files) {
    FastqReader* reader {new FastqReader(file.c_str())};
    readers.push_back(reader);
  }

  ReadVector_ readVector {readFastq_(readers)};
  while (not readVector.eof) {
    co_yield readVector.reads;
    readVector.free();
    readVector = readFastq_(readers);
  }
  // TODO: Extra readVector.free() here when files are not of equal length?

  for (FastqReader* const reader: readers) {
    delete reader;
  }
}

vector<char> getNucleotides(
    vector<Read*> const& reads, vector<size_t> const ntToTake,
    size_t const headerUMISize) {
  vector<char> nucleotides;

  // Pull the UMI from the header of the first read.
  if (headerUMISize > 0) {
    // Get the UMI, and cut/pad it to headerUMISize.
    string headerUMI {makeStringSize_(
      extractUMI(reads.front()), headerUMISize, 'N')};
    for (size_t i {0}; i < headerUMISize; i++) {
      nucleotides.push_back(headerUMI[i]);
    }
  }

  for (size_t i {0}; i < reads.size(); i++) {
    Read* read {reads[i]};
    size_t length {ntToTake[i]};

    // Padd the reads with N if it is too short.
    string sequence {makeStringSize_(*read->mSeq, length, 'N')};

    // Add length nucleotides from Read.
    for (size_t pos {0}; pos < length; pos++) {
      nucleotides.push_back(sequence[pos]);
    }
  }
  return nucleotides;
}

Word makeWord(
    vector<Read*> const& reads, vector<size_t> const ntToTake,
    size_t const headerUMISize) {
  Word word;
  vector<char> nucleotides {getNucleotides(reads, ntToTake, headerUMISize)};
  for (char const& nucleotide: nucleotides) {
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

void printWord(vector<uint8_t> const& word) {
  for (uint8_t const& letter: word) {
    cout << ' ' << (int)letter;
  }
  cout << '\n';
}

string addDir(char const filename[], string const dir) {
  return dir + '/' + filename;
}

string makeFileName(
    string const filename, string const dir, string const suffix) {
  string name {basename(filename)};
  size_t pos {name.find('.')};
  string suff {name.substr(0, pos) + '_' + suffix +
    name.substr(pos, string::npos)};
  return addDir(suff.c_str(), dir);
}

vector<string> makeFileNames(
    vector<string> const files, string const dir, string const suffix) {
  vector<string> fileNames;
  for (string const& name: files) {
    fileNames.push_back(makeFileName(name, dir, suffix));
  }
  return fileNames;
}

string extractLastField(string const str, char const sep) {
  size_t last {str.find_last_of(sep)};

  if (last != string::npos) {
    return str.substr(last + 1);
  }
  return "";
}

bool validUMI(string const umi) {
  // An empty UMI is not valid.
  if (umi.empty()) {
    return false;
  }

  // Only ATCG is valid in a UMI.
  for (char const c: umi) {
    if (nuc.find(c) == nuc.end()) {
      return false;
    }
  }
  return true;
}

string extractUMI(Read* const read) {
  return extractUMI_(*read->mName);
}

vector<size_t> ntFromFile(size_t const files, size_t const length) {
  vector<size_t> v{};
  size_t div {length / files};
  // All items are set to div, except the last one.
  for (size_t i {0}; i < files - 1; i++) {
    v.push_back(div);
  }
  // Remainder is added to the last item.
  v.push_back(div + length % files);
  return v;
}
