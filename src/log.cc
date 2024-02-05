#include "log.h"


time_t startMessage(ofstream& log, char const message[]) {
  log << message << "... ";
  log.flush();

  return time(nullptr);
}

void endMessage(ofstream& log, time_t const start) {
  time_t seconds {static_cast<time_t>(difftime(time(nullptr), start))};
  log << "done. (" << seconds / 60 << 'm' << seconds % 60 << "s)\n";
  log.flush();
}
