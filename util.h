/**
 * @file util.h
 * @author László Veréb
 * @date 2010.03.27
 */

#ifndef UTIL_H
#define UTIL_H

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#ifndef NDEBUG
#include <string.h>
#endif

#define CALLER_LENGTHS 100

extern char caller_Function[CALLER_LENGTHS];
extern char caller_Function_Pretty[CALLER_LENGTHS];
extern char caller_File[CALLER_LENGTHS];
extern short caller_Line;

//typedef double * const dpc;

/**	logikai változótípus definiálása ismeretlen értékkel is
 */
typedef enum {
	false = 0, true = 1
} bool;

/**
 *		The function is a secured variant of the built-in fopen.
 * @param[in]	name	: name of the file
 * @param[in]	mode	: with wich mode to open it
 * @return	pointer to the opened file
 */
FILE * sfopen(char *name, char * mode);

/**
 *		The function is a secured read-only variant of the built-in fopen.
 * @param[in]	name	: name of the file
 * @return	pointer to the opened file
 */
FILE * sfopen_read(char *name);

/**
 *		The function is a secured write-only variant of the built-in fopen.
 * @param[in]	name	: name of the file
 * @return	pointer to the opened file
 */
FILE * sfopen_write(char *name);

#endif // UTIL_H
