/**
 * @file main_test.c
 *
 * @date Jul 30, 2011
 * @author vereb
 * @brief
 */

#include <time.h>
#include <stdlib.h>
#include "util_math.h"
#include "binary_system.h"

int main(int argc, char *argv[]) {
	argc = argc;
	argv = argv;
	srand(time(NULL));
#ifdef TEST
	areUtilMathFunctionsOK();
	areIOFunctionsGood();
	areBinarySystemMassFunctionsGood();
	areBinarySystemSpinFunctionsGood();
#endif // TEST
	puts("\nOK");
	return 0;
}
