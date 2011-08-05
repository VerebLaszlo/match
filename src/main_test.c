/**
 * @file main_test.c
 *
 * @date Jul 30, 2011
 * @author vereb
 * @brief
 */

#include <math.h>
#include <time.h>
#include <stdlib.h>
#include "util_IO.h"
#include "util_math.h"
#include "binary_system.h"

int main(int argc, char *argv[]) {
	argc = argc;
	argv = argv;
	srand(time(NULL));
	areUtilMathFunctionsOK();
	areIOFunctionsGood();
	areBinarySystemMassFunctionsGood();
	areBinarySystemSpinFunctionsGood();
	/*OutputFormat format;
	setOutputFormat(&format, 5, 5, '%', true, "plot", 0);
	binarySystem system, limits[2];
	for (ushort i = 0; i < 2; i++) {
		double j = i + 1;
		limits[i].mass.mass[0] = j * 30.0;
		limits[i].mass.mass[1] = j * 3.0;
		limits[i].spin[0].magnitude = j * 0.5;
		limits[i].spin[0].azimuth[FIXED] = j * M_PI / 3.0;
		limits[i].spin[0].inclination[FIXED] = j * M_PI / 4.0;
	}
	generateBinarySystemParameters(&system, limits, FROM_M1M2, FROM_FIXED_ANGLES);
	printBinarySystemParameters(stdout, &system, &format);*/
	puts("\nOK");
	return 0;
}
