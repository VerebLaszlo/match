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
