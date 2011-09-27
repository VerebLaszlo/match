/**
 * @file main_test.c
 *
 * @date Jul 30, 2011
 * @author vereb
 * @brief
 */

#include "program_functions.h"
#include <string.h>
#include <stdlib.h>

static void printHelp(void) {
	puts("h: prints this help");
	puts("p: save the generated waveform for plotting");
	puts("m: calculate the match");
	puts("s: the file containing the parameters");
	puts("i: the file initialize the program");
}

typedef struct {
	char parameterFile[FILENAME_MAX];
	char programFile[FILENAME_MAX];
	bool plot;
	bool calculateMatch;
} Options;

static void interpretOptions(Options *option, int argc, char *argv[]) {
	option->plot = false;
	option->calculateMatch = false;
	while (argc) {
		if ((*argv)[0] == '-') {
			if ((*argv)[1] == 'h') {
				printHelp();
				exit(EXIT_SUCCESS);
			} else if ((*argv)[1] == 'p') {
				option->plot = true;
			} else if ((*argv)[1] == 'm') {
				option->calculateMatch = true;
			} else if ((*argv)[1] == 's') {
				argc--;
				argv++;
				strcpy(option->parameterFile, *argv);
			} else if ((*argv)[1] == 'i') {
				argc--;
				argv++;
				strcpy(option->programFile, *argv);
			}
		}
		argc--;
		argv++;
	}
}

int main(int argc, char *argv[]) {
#ifdef TEST
	//testingFunctions();
#endif // TEST
	Options option;
	if (argc == 1) {
		printHelp();
		exit(EXIT_SUCCESS);
	}
	argc--;
	argv++;
	interpretOptions(&option, argc, argv);
	run(option.programFile, option.parameterFile, option.plot, option.calculateMatch);
	puts("\nOK");
	return 0;
}
