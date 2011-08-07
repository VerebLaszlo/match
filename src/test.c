/**
 * @file test.c
 * @author László Veréb
 * @date 2011.08.07.
 * @brief Macros for testing functions.
 */

#include "test.h"

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
