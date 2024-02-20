#pragma once

#include <string>
#include <vector>

#include "../lib/trie/lib/CPP20Coroutines/include/generator.hpp"
#include "../lib/fastp/src/fastqreader.h"

using std::string;
using std::vector;

/*! Word. */
struct Word {
  vector<uint8_t> data {};
  bool filtered {false};
};


/*! Loop over all reads in multiple FastQ files.
 *
 * \param files FastQ file names.
 *
 * \return All reads.
 */
generator<vector<Read*>> readFiles(vector<string> const);

/*! Extract `wordLength` nucleotides from `reads`. If the first file has a UMI
 * in the header, this will get preference.
 */
vector<char> getNucleotides(
  vector<Read*> const&, vector<size_t> const, size_t const);

/*! Select a total of `wordLength` nucleotides from every read in `reads` to
 * create a word.
 *
 * \param reads Reads.
 * \param wordLength Read selection length.
 *
 * \return Word.
 */
Word makeWord(vector<Read*> const&, vector<size_t> const, size_t const);

/*! Print a word.
 *
 * \param word Word.
 */
void printWord(vector<uint8_t> const&);

/*!
 */
string addDir(char const[], string const);

/*!
 */
string makeFileName(string const, string const, string const);

/*!
 */
vector<string> makeFileNames(vector<string> const, string const, string const);

/*! Extract the last field from string.
 *
 * \param str String to extract field from.
 * \param sep Separator between the fields.
 */
string extractLastField(string const, char const);

/*! Determine of UMI is a valid UMI. It must be non-emtpy and only contain
 * characters from ATCG.
 *
 * \param umi The UMI to check.
 */
bool validUMI(string const);

/*! Extract UMI from a read.
 *
 * \param read Read.
 */
string extractUMI(Read* const);

/*! Divide `length` nucleotides over `files`, with the remainder used on the
 * last file.
 *
 * \param files Number of files.
 * \param lengt_h Total number of nucleotides to divide.
 */
vector<size_t> ntFromFile(size_t const, size_t const);
