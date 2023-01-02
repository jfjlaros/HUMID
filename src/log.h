#pragma once

#include <fstream>

using std::ofstream;


/*!
 * Write a task start message to a log.
 *
 * \param log Log file.
 * \param message Message.
 *
 * \return Task start time.
 */
time_t startMessage(ofstream&, char const[]);

/*!
 * Write a task end message to a log.
 *
 * \param log Log file.
 * \param start Task start time.
 */
void endMessage(ofstream&, time_t const);

/*!
 * \param value
 */
char const* boolToText(bool const);
