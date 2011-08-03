/**
 * @file util.c
 * @author László Veréb
 * @date 2011.07.19.
 * @brief Contains useful functions.
 */

#include "util.h"

void neg(bool *var) {
	*var = !*var;
}

#ifdef TEST
_sourceLocation sourceLocation;
char previous_function[FILENAME_MAX];
char function_file[FILENAME_MAX];
ushort function_line;
char normal[8] = "\e[0m";
char err[8] = "\e[0;31m";
char errBold[8] = "\e[1;31m";
char ok[8] = "\e[0;36m";
char okBold[8] = "\e[1;36m";
#endif
