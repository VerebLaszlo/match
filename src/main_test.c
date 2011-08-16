/**
 * @file main_test.c
 *
 * @date Jul 30, 2011
 * @author vereb
 * @brief
 */

#include <time.h>
#include <stdlib.h>
#include <math.h>
#include "util_math.h"
#include "binary_system.h"
#include "detector.h"

int main(int argc, char *argv[]) {
	printf("%d: %s\n", argc, argv[0]);
	srand(time(NULL));
#ifdef TEST
	areUtilMathFunctionsOK();
	areIOFunctionsGood();
	areBinarySystemMassFunctionsGood();
	areBinarySystemSpinFunctionsGood();
	areBinarySystemFunctionsGood();
	areDetectorFunctionsGood();
#else
	bool boolean = true;
	printf("%d - ", boolean);
	neg(&boolean);
	printf("%d\n", boolean);
	printf("%lg, %lg, %lg\n", randomBetweenZeroAndOne(), randomBetweenZeroAnd(100.0),
		randomBetween(-100.0, 100.0));
	printf("%lg, %lg, %lg\n", sinGood(M_PI_4), cosGood(M_PI_4 + M_PI_2), tanGood(M_PI_4));
	printf("%lg, %lg, %lg\n", radianFromDegree(90.0), radianFromTurn(0.5), turnFromDegree(90.0));
	printf("%lg, %lg, %lg\n", degreeFromRadian(M_PI_4), degreeFromTurn(0.75), turnFromRadian(M_PI));
	double s = square(2.0);
	double c = cube(3.0);
	printf("%lg, %lg, %lg\n", normaliseRadians(5.0 * M_PI_2), s, c);
	printf("%d\n", isNear(1.0, 2.0, 1.1));
	OutputFormat format;
	setOutputFormat(&format, 5, 11, '%', true, "x", 1);
	short length = 2 * format.widthWithSeparator;
	char formatString[length];
	setFormat(formatString, 2, &format);
	char formatEndString[length];
	setFormatEnd(formatEndString, 2, &format);
	printf("%s %s %s", format.oneNumber, formatString, formatEndString);
	massParameters mass, massLimits[2];
	massLimits[MAX].mass[0] = 10.0 * (massLimits[MIN].mass[0] = 2.0);
	massLimits[MAX].mass[1] = 10.0 * (massLimits[MIN].mass[1] = 3.0);
	massLimits[MAX].totalMass = 10.0 * (massLimits[MIN].totalMass = 5.0);
	massLimits[MAX].mu = 10.0 * (massLimits[MIN].mu = 1.0);
	massLimits[MAX].eta = 2.5 * (massLimits[MIN].eta = 0.1);
	massLimits[MAX].chirpMass = -1.0 * (massLimits[MIN].chirpMass = -10.0);
	massLimits[MAX].nu = -1.0 * (massLimits[MIN].nu = -10.0);
	massLimits[MAX].m1_m2 = -1.0 * (massLimits[MIN].m1_m2 = -10.0);
	generateMass(&mass, massLimits, GEN_M1M2);
	printMassParameters(stdout, &mass, &format);
	spinParameters spin[2], spinLimits[2];
	spinLimits[MAX].magnitude = 10.0 * (spinLimits[MIN].magnitude = 2.0);
	for (short coordinate = 0; coordinate < COORDINATE_CONVENTIONS; coordinate++) {
		for (short dim = 0; dim < DIMENSION; dim++) {
			spinLimits[MAX].component[coordinate][dim] = -1.0
				* (spinLimits[MIN].component[coordinate][dim] = -10.0);
			spinLimits[MAX].unity[coordinate][dim] = -1.0
				* (spinLimits[MIN].unity[coordinate][dim] = -1.0);
		}
		spinLimits[MAX].azimuth[coordinate] = 600.0 * (spinLimits[MIN].azimuth[coordinate] = 0.01);
		spinLimits[MAX].inclination[coordinate] = 300.0 * (spinLimits[MIN].inclination[coordinate] =
			0.01);
		spinLimits[MAX].elevation[coordinate] = -1.0
			* (spinLimits[MIN].elevation[coordinate] = -1.0);
	}
	generateSpin(&spin[0], spinLimits, 1.0, GEN_FIXED_XYZ);
	printSpinParameters(stdout, &spin[0], &format);
	printSpinParameters(stdout, &spin[1], &format);
#endif // TEST
	puts("\nOK");
	return 0;
}
