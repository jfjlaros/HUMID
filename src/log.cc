#include "log.h"


/*!
 * Write a task start message to a log.
 *
 * \param log Log file.
 * \param message Message.
 *
 * \return Task start time.
 */
time_t startMessage(ofstream& log, char const message[]) {
  log << message << "... ";
  log.flush();

  return time(nullptr);
}

/*!
 * Write a task end message to a log.
 *
 * \param log Log file.
 * \param start Task start time.
 */
void endMessage(ofstream& log, time_t const start) {
  time_t seconds {static_cast<time_t>(difftime(time(nullptr), start))};
  log << "done. (" << seconds / 60 << 'm' << seconds % 60 << "s)\n";
  log.flush();
}

char const* boolToText(bool const value) {
  if (value) {
    return "true";
  }
  return "false";
}
