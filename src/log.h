#pragma once

#include <fstream>

using std::ofstream;

time_t startMessage(ofstream&, char const[]);
void endMessage(ofstream&, time_t const);

/*!
 * \param value
 */
char const* boolToText(bool const);
