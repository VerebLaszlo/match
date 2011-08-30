/**
 * @file binary_system.c
 *
 * @date Jul 20, 2011
 * @author vereb
 * @brief Binary system specific
 */

#include "test.h"
#include <assert.h>
#include "util_math.h"
#include "binary_system.h"

/// @name Generation functions
///@{

void generateBinarySystemParameters(BinarySystem *system, BinarySystem limits[],
	generationMode genMass, generationMode genSpin) {
	BACKUP_DEFINITION_LINE(); //
	assert(system);
	assert(limits);
	system->inclination = randomBetween(limits[MIN].inclination, limits[MAX].inclination);
	system->distance = randomBetween(limits[MIN].distance, limits[MAX].distance);
	for (short i = 0; i < NUMBER_OF_BLACKHOLES; i++) {
		system->flatness[i] = randomBetween(limits[MIN].flatness[i], limits[MAX].flatness[i]);
	}
	massParameters mass[2] = { limits[MIN].mass, limits[MAX].mass };
	generateMass(&system->mass, mass, genMass);
	spinParameters spin[2];
	for (short i = 0; i < NUMBER_OF_BLACKHOLES; i++) {
		spin[MIN] = limits[MIN].spin[i];
		spin[MAX] = limits[MAX].spin[i];
		generateSpin(&system->spin[i], spin, system->inclination, genSpin);
	} //
	SAVE_FUNCTION_FOR_TESTING();
}

///@}
/// @name Printing functions
///@{

void printBinarySystemParameters(FILE *file, BinarySystem *system, OutputFormat *format) {
	BACKUP_DEFINITION_LINE();
	printMassParameters(file, &system->mass, format);
	printSpinParameters(file, &system->spin[0], format);
	printSpinParameters(file, &system->spin[1], format);
	ushort number = 3;
	ushort length = number * format->widthWithSeparator;
	char formatString[length];
	setFormat(formatString, number, format);
	fprintf(file, formatString, system->flatness[0], system->flatness[1], system->inclination);
	setFormatEnd(formatString, number, format);
	fprintf(file, formatString, system->distance, system->coalescencePhase,
		system->coalescenceTime);
	SAVE_FUNCTION_FOR_TESTING();
}

///@}

//#ifdef TEST
#include <math.h>
/// @name Testing functions
///@{

static bool isOK_generateBinarySystemParameters(void) {
	bool isOK = true;
	if (!areUtilMathFunctionsOK()) {
		isOK = false;
	} else if (!areBinarySystemMassFunctionsGood()) {
		isOK = false;
	} else if (!areBinarySystemSpinFunctionsGood()) {
		isOK = false;
	}
	BinarySystem system, limits[2];
	SAVE_FUNCTION_CALLER();
	limits[MAX].mass.mass[0] = 100.0 * (limits[MIN].mass.mass[0] = 1.0);
	limits[MAX].mass.mass[1] = 100.0 * (limits[MIN].mass.mass[1] = 1.0);
	limits[MAX].mass.totalMass = 100.0 * (limits[MIN].mass.totalMass = 2.0);
	limits[MAX].mass.eta = 100.0 + (limits[MIN].mass.eta = 0.0);
	limits[MAX].mass.chirpMass = 100.0 + (limits[MIN].mass.chirpMass = 0.0);
	limits[MAX].mass.mu = 100.0 + (limits[MIN].mass.mu = 0.0);
	limits[MAX].mass.nu = 100.0 + (limits[MIN].mass.nu = 0.0);
	limits[MAX].mass.m1_m2 = 100.0 + (limits[MIN].mass.m1_m2 = 0.0);
	for (ushort dim = 0; dim < DIMENSION; dim++) {
		limits[MAX].spin[0].component[FIXED][dim] = 100.0
			* (limits[MIN].spin[0].component[FIXED][dim] = 1.0);
		limits[MAX].spin[1].component[FIXED][dim] = 100.0
			* (limits[MIN].spin[1].component[FIXED][dim] = 1.0);
	}
	for (ushort i = 0; i < 2; i++) {
		limits[MAX].spin[i].magnitude = 100.0 * (limits[MIN].spin[i].magnitude = 1.0);
		limits[MAX].spin[i].azimuth[FIXED] = 10000.0 * (limits[MIN].spin[i].azimuth[FIXED] =
			0.000314);
		limits[MAX].spin[i].azimuth[PRECESSING] = 10000.0
			* (limits[MIN].spin[i].azimuth[PRECESSING] = 0.000314);
		limits[MAX].spin[i].inclination[FIXED] = 10000.0 * (limits[MIN].spin[i].inclination[FIXED] =
			0.000314);
		limits[MAX].spin[i].inclination[PRECESSING] = 10000.0
			* (limits[MIN].spin[i].inclination[PRECESSING] = 0.000314);
		limits[MAX].spin[i].elevation[FIXED] = M_PI + (limits[MIN].spin[i].elevation[FIXED] =
			-M_PI_2);
		limits[MAX].spin[i].elevation[PRECESSING] = M_PI
			+ (limits[MIN].spin[i].elevation[PRECESSING] = -M_PI_4);
		limits[MAX].spin[i].component[FIXED][X] = 10.0 * (limits[MIN].spin[i].component[FIXED][X] =
			1.0);
		limits[MAX].spin[i].component[FIXED][Y] = 10.0 * (limits[MIN].spin[i].component[FIXED][Y] =
			1.0);
		limits[MAX].spin[i].component[FIXED][Z] = 10.0 * (limits[MIN].spin[i].component[FIXED][Z] =
			1.0);
		limits[MAX].spin[i].component[PRECESSING][X] = 40.0
			+ (limits[MIN].spin[i].component[PRECESSING][X] = -20.0);
		limits[MAX].spin[i].component[PRECESSING][Y] = 40.0
			+ (limits[MIN].spin[i].component[PRECESSING][Y] = -20.0);
		limits[MAX].spin[i].component[PRECESSING][Z] = 40.0
			+ (limits[MIN].spin[i].component[PRECESSING][Z] = -20.0);
	}
	generateBinarySystemParameters(&system, limits, GEN_M1M2, GEN_FIXED_XYZ);
	if (!isOK) {
		return isOK;
	}PRINT_OK();
	return isOK;
}

bool areBinarySystemFunctionsGood(void) {
	bool isOK = true;
	if (!isOK_generateBinarySystemParameters()) {
		isOK = false;
	}
	if (isOK) {
		PRINT_OK_FILE();
	} else {
		PRINT_ERROR_FILE();
	}
	return isOK;
}

///@}
//#endif	// TEST
