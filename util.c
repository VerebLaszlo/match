/*
 * @file util.c
 * @author László Veréb
 * @date 2010.01.11.
 */

#include "util.h"

extern char * program_invocation_short_name;
extern char * program_invocation_name;

FILE * safely_Open_File(char *file_Name, char *mode) {
	assert(strcmp(file_Name, ""));
	assert(strcmp(mode, ""));
	FILE *stream;
	errno = 0;
	stream = fopen (file_Name, mode);
	if (stream == NULL) {
		fprintf (stderr, "%s: Couldn't open file %s; %s\n", program_invocation_short_name, file_Name, strerror (errno));
		exit (EXIT_FAILURE);
	} else return stream;
}

FILE * safely_Open_File_For_Reading(char *file_Name) {
	assert(strcmp(file_Name, ""));
	FILE *stream;
	errno = 0;
	stream = fopen (file_Name, "r");
	if (stream == NULL) {
		fprintf (stderr, "%s: Couldn't open file for reading %s; %s\n", program_invocation_short_name, file_Name, strerror (errno));
		exit (EXIT_FAILURE);
	} else return stream;
}

FILE * safely_Open_File_For_Writing(char *file_Name) {
	assert(strcmp(file_Name, ""));
	FILE *stream;
	errno = 0;
	stream = fopen (file_Name, "w");
	if (stream == NULL) {
		fprintf (stderr, "%s: Couldn't open file for writing %s; %s\n", program_invocation_short_name, file_Name, strerror (errno));
		exit (EXIT_FAILURE);
	} else return stream;
}
