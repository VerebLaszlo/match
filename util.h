/**
 * @file util.h
 * @author László Veréb
 * @date 2010.03.27
 */

#ifndef UTIL_H_
#define UTIL_H_

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <math.h>

typedef double * const dpc;

//	logikai változótípus definiálása ismeretlen értékkel is
typedef enum {
	null = -1, false = 0, true = 1
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

//	egyéb

/**
 *		The function generates random number in the [0,1) intervall.
 * @return	a random number
 */
double rand1(void);

/**
 *		The function rounds the number to the closest power of two.
 * @param[in]	num	: number to round
 * @return	rounded value
 */
long round_po2(double num);

#endif /* UTIL_H_ */
