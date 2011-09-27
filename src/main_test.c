/**
 * @file main_test.c
 *
 * @date Jul 30, 2011
 * @author vereb
 * @brief
 */

#include <time.h>
#include <stdlib.h>
#include "program_functions.h"

int main(int argc, char *argv[]) {
	printf("%d: %s\n", argc, argv[0]);
#ifdef TEST
	srand(86);
	areUtilMathFunctionsOK();
	areIOFunctionsGood();
	areBinarySystemMassFunctionsGood();
	areBinarySystemSpinFunctionsGood();
	areBinarySystemFunctionsGood();
	areDetectorFunctionsGood();
#endif // TEST
	run(argv[1], argv[2], true, false);
	puts("\nOK");
	return 0;
}
