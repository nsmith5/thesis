#include <signal.h>

extern volatile sig_atomic_t time_to_leave;

void init (int argc, char **argv);

void finalize (void);
