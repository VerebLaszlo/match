/*
 * @file util.c
 * @author László Veréb
 * @date 2010.01.11.
 */

#include "util.h"

extern char * program_invocation_short_name;
extern char * program_invocation_name;

FILE * sfopen(char *name, char *mode) {
	FILE *stream;
	errno = 0;
	stream = fopen (name, mode);
	if (stream == NULL) {
		fprintf (stderr, "%s: Couldn't open file %s; %s\n", program_invocation_short_name, name, strerror (errno));
		exit (EXIT_FAILURE);
	} else return stream;
}

FILE * sfopen_read(char *name) {
	FILE *stream;
	errno = 0;
	stream = fopen (name, "r");
	if (stream == NULL) {
		fprintf (stderr, "%s: Couldn't open file for reading %s; %s\n", program_invocation_short_name, name, strerror (errno));
		exit (EXIT_FAILURE);
	} else return stream;
}

FILE * sfopen_write(char *name) {
	FILE *stream;
	errno = 0;
	stream = fopen (name, "w");
	if (stream == NULL) {
		fprintf (stderr, "%s: Couldn't open file for writing %s; %s\n", program_invocation_short_name, name, strerror (errno));
		exit (EXIT_FAILURE);
	} else return stream;
}

//	egyéb

double rand1(void) {
	return ((double)rand()) / ((double)RAND_MAX);
}

long round_po2(double num) {
	register double temp = log(num) / M_LN2;
#ifdef __USE_ISOC99

	return (long) exp2(ceil(temp));
#endif
	//return (long) pow(2., ceil(temp));
	return (long) exp(ceil(temp) * M_LN2);
}
