#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
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
  nucleotides.reserve(wordLength);

  // Pull the UMI from the header of the first read
  string headerUMI = extractUMI(reads.front());
  for (size_t i = 0; i < wordLength and i < headerUMI.size(); i++) {
      nucleotides.push_back(headerUMI[i]);
  }

  // The length we still have available from wordLength after extracting the
  // UMI from the header
  size_t length = wordLength - nucleotides.size();

  // The number of nucleotides to take from each file
  vector<size_t> ntToTake = ntFromFile(reads.size(), length);

  for (size_t i = 0; i < reads.size(); i++) {
    Read* read = reads[i];
    size_t length = ntToTake[i];

    // Throw if the read is too short
    if ((*read->mSeq).size() < length) {
      std::ostringstream msg;
      msg << "Attempted to read "
          << length
          << " nucleotides from "
          << *read->mName
          << " (length="
          << (*read->mSeq).size()
          << ").";
      throw std::out_of_range( msg.str() );
    }

    // Add length nucleotides from Read
    for (size_t pos = 0; pos < length; pos++) {
      char nucleotide = (*read->mSeq)[pos];
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

  // The UMI must be before the first space
  string substr = header.substr(0, first_space);

  // If we detect a UMI with a _ separator, we return that UMI
  string UMI = extractLastField(substr, '_');

  // Check if the UMI only contains characters from 'ATCGN'. If we find any
  // other character, the 'UMI' is not actually a UMI.
  if (validUMI(UMI)) {
    return UMI;
  }

  // Otherwise, we check to see if we find a BCL Convert style UMI
  UMI = extractLastField(substr, ':');
  if (validUMI(UMI)) {
    return UMI;
  }
  else {
    return "";
  }
}

/*
 * Extract the last field from string
 *
 * \param string String to extract field from
 * \param sep Separator between the fields
 */
string extractLastField(string string, char sep) {
  size_t last = string.find_last_of(sep);

  if (last != string::npos) {
    return string.substr(last + 1);
  }
  else {
    return "";
  }
}


/*!
 * Determine of UMI is a valid UMI. It must be non-emtpy and only contain
 * characters from ATCGN
 *
 * \param UMI The UMI to check
 */
bool validUMI(string UMI) {
  // An empty UMI is not valid
  if (UMI.empty()) {
    return false;
  }

  // Only ATCGN is valid in a UMI
  for (char c: UMI) {
    if (nnuc.find(c) == nnuc.end()) {
      return false;
    }
  }
  return true;
}

/*!
 * Extract UMI from a read
 *
 * \param read Read
 */
string extractUMI(Read* read) {
  return _extractUMI(*read->mName);
}

/*!
 * Divide `length` nucleotides over `files`, with the remainder used on the
 * last file.
 *
 * \param files Number of files.
 * \param lengt_h Total number of nucleotides to divide.
 */
vector<size_t> ntFromFile(size_t files, size_t length) {
  vector<size_t> v{};
  size_t div = length / files;
  size_t remainder = length % files;
  // All items are set to div, except the last one
  for (size_t i = 0; i < files - 1; i++) {
    v.push_back(div);
  }
  // Remainder is added to the last item
  v.push_back(div+remainder);
  return v;
}
