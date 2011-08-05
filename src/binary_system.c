/**
 * @file binary_system.c
 *
 * @date Jul 20, 2011
 * @author vereb
 * @brief Binary system specific
 */

#include "binary_system.h"

/// @name Helper functions
///@{

/** Returns the magnitude of the spins.
 * @param[in,out] spin	: spin components
 */
static void magnitudeOfSpin(spinParameters *spin) {
	BACKUP_DEFINITION_LINE(); //
	assert(spin);
	spin->magnitude = 0.0;
	for (short j = 0; j < DIMENSION; j++) {
		spin->magnitude += square(spin->component[FIXED][j]);
	}
	spin->magnitude = sqrt(spin->magnitude);
	SAVE_FUNCTION_FOR_TESTING();
}

/**	Returns true, if the mass parameters are between their limits.
 * @param[in] spin		: spin parameters to examine
 * @param[in] limits	: limits of the spin parameters
 * @return true or false
 */
static bool isSpinBetweenLimits(spinParameters *spin, spinParameters limits[]) {
	BACKUP_DEFINITION_LINE(); //
	assert(spin);
	assert(limits);
	if (limits[MIN].magnitude > spin->magnitude || spin->magnitude > limits[MAX].magnitude) {
		return false;
	}
	for (short coordinate = 0; coordinate < COORDINATE_CONVENTIONS; coordinate++) {
		if (limits[MIN].azimuth[coordinate] > spin->azimuth[coordinate]
				|| spin->azimuth[coordinate] > limits[MAX].azimuth[coordinate]) {
			return false;
		}
		if (limits[MIN].inclination[coordinate] > spin->inclination[coordinate]
				|| spin->inclination[coordinate] > limits[MAX].inclination[coordinate]) {
			return false;
		}
		if (limits[MIN].elevation[coordinate] > spin->elevation[coordinate]
				|| spin->elevation[coordinate] > limits[MAX].elevation[coordinate]) {
			return false;
		}
		for (short dimension = 0; dimension < DIMENSION; dimension++) {
			if (limits[MIN].component[coordinate][dimension]
					> spin->component[coordinate][dimension]
					|| spin->component[coordinate][dimension]
							> limits[MAX].component[coordinate][dimension]) {
				return false;
			}
		}
	}
	SAVE_FUNCTION_FOR_TESTING();
	return true;
}

///@}
/// @name Conversion functions
///@{

/**	Converts spin components to angles in the specified frame.
 * @param[in,out]	spin	: spin parameters
 * @param[in]		index	: frame
 */
static void convertSpinFromXyzToAngles(spinParameters *spin, const ushort index) {
	BACKUP_DEFINITION_LINE();
	spin->inclination[index] = acos(spin->component[index][Z] / spin->magnitude);
	spin->azimuth[index] = -acos(
			spin->component[index][X] / spin->magnitude / sin(spin->inclination[index]));
	if (spin->component[index][Y] < 0.0) {
		spin->azimuth[index] *= -1.0;
	}
	if (!isfinite(spin->azimuth[index])) {
		spin->azimuth[index] = 0.0;
	}
	SAVE_FUNCTION_FOR_TESTING();
}

/**	Converts spin angles to components in the specified frame.
 * @param[in,out]	spin	: spin parameters
 * @param[in]		index	: frame
 */
static void convertSpinFromAnglesToXzy(spinParameters *spin, const ushort index) {
	BACKUP_DEFINITION_LINE();
	double cosAzimuth, cosInclination;
	cosAzimuth = spin->azimuth[index] == M_PI_2 ? 0.0 : cos(spin->azimuth[index]);
	cosInclination = spin->inclination[index] == M_PI_2 ? 0.0 : cos(spin->inclination[index]);
	spin->component[index][X] = spin->magnitude * sin(spin->inclination[index]) * cosAzimuth;
	spin->component[index][Y] = spin->magnitude * sin(spin->inclination[index])
			* sin(spin->azimuth[index]);
	spin->component[index][Z] = spin->magnitude * cosInclination;
	SAVE_FUNCTION_FOR_TESTING();
}

/**	Converts spin parameters from the fixed frame to the precessing frame.
 * @param[in,out]	spin		: spin parameters
 * @param[in]		inclination	: inclination of the precessing frame
 */
static void convertSpinFromFixedFrame(spinParameters *spin, const double inclination) {
	BACKUP_DEFINITION_LINE(); //
	assert(spin);
	double theta[2] = { +0.0, inclination };
	spin->component[PRECESSING][X] = spin->component[FIXED][X] * cos(theta[0]) * cos(theta[1])
			+ spin->component[FIXED][Y] * sin(theta[0]) * cos(theta[1])
			- spin->component[FIXED][Z] * sin(theta[1]);
	spin->component[PRECESSING][Y] = -spin->component[FIXED][X] * sin(theta[0])
			+ spin->component[FIXED][Y] * cos(theta[0]);
	spin->component[PRECESSING][Z] = spin->component[FIXED][X] * cos(theta[0]) * sin(theta[1])
			+ spin->component[FIXED][Y] * sin(theta[0]) * sin(theta[1])
			+ spin->component[FIXED][Z] * cos(theta[1]);
	SAVE_FUNCTION_FOR_TESTING();
}

/**	Converts spin parameters from the precessing frame to the fixed frame.
 * @param[in,out]	spin		: spin parameters
 * @param[in]		inclination	: inclination of the precessing frame
 */
static void convertSpinFromPrecessionFrame(spinParameters *spin, const double inclination) {
	BACKUP_DEFINITION_LINE(); //
	assert(spin);
	double theta[2] = { -0.0, inclination };
	spin->component[FIXED][X] = spin->magnitude
			* (+spin->component[PRECESSING][X] * cos(theta[0]) * cos(theta[1])
					+ spin->component[PRECESSING][Y] * sin(theta[0])
					- spin->component[PRECESSING][Z] * cos(theta[0]) * sin(theta[1]));
	spin->component[FIXED][Y] = spin->magnitude
			* (-spin->component[PRECESSING][X] * sin(theta[0]) * sin(theta[1])
					+ spin->component[PRECESSING][Y] * cos(theta[0])
					+ spin->component[PRECESSING][Z] * sin(theta[0]) * sin(theta[1]));
	spin->component[FIXED][Z] = spin->magnitude
			* (+spin->component[PRECESSING][X] * sin(theta[1])
					+ spin->component[PRECESSING][Y] * cos(theta[1]));
	SAVE_FUNCTION_FOR_TESTING();
}

/**	Converts spin parameters according the conversion parameter.
 * @param[in,out]	spin		: spin parameters
 * @param[in]		inclination	: inclination fo the precessing frame
 * @param[in]		convert		: conversion mode
 */
static void convertSpin(spinParameters spin[NUMBER_OF_BLACKHOLES], const double inclination,
		conversionMode convert) {
	BACKUP_DEFINITION_LINE(); //
	assert(spin);
	switch (convert) {
	case FROM_FIXED_XYZ:
		magnitudeOfSpin(spin);
		convertSpinFromFixedFrame(spin, inclination);
		convertSpinFromXyzToAngles(spin, FIXED);
		convertSpinFromXyzToAngles(spin, PRECESSING);
		break;
	case FROM_FIXED_ANGLES:
		convertSpinFromAnglesToXzy(spin, FIXED);
		convertSpinFromFixedFrame(spin, inclination);
		convertSpinFromXyzToAngles(spin, PRECESSING);
		break;
	case FROM_PRECESSION_XZY:
		magnitudeOfSpin(spin);
		convertSpinFromXyzToAngles(spin, FIXED);
		convertSpinFromPrecessionFrame(spin, inclination);
		convertSpinFromXyzToAngles(spin, PRECESSING);
		break;
	case FROM_PRECESSION_ANGLES:
		convertSpinFromAnglesToXzy(spin, PRECESSING);
		convertSpinFromPrecessionFrame(spin, inclination);
		convertSpinFromXyzToAngles(spin, FIXED);
		break;
	default:
		break;
	}
	SAVE_FUNCTION_FOR_TESTING();
}

///@}
/// @name Generation functions
///@{

/** Generates the spin parameters according the generation mode.
 * @param[out]	spin		: generated spin parameters
 * @param[in]	limits		: limits of the spin parameters
 * @param[in]	inclination	: inclination of the precessing frame
 * @param[in]	mode		: generation mode
 */
static void generateSpin(spinParameters *spin, spinParameters limits[], double inclination,
		conversionMode mode) {
	BACKUP_DEFINITION_LINE(); //
	assert(spin);
	assert(limits);
	switch (mode) {
	case GEN_FIXED_XYZ:
		for (short i = 0; i < NUMBER_OF_BLACKHOLES; i++) {
			do {
				for (short j = 0; j < DIMENSION; j++) {
					spin[i].component[FIXED][j] = randomBetween(limits[MIN].component[FIXED][j],
							limits[MAX].component[FIXED][j]);
				}
				convertSpin(&spin[i], inclination, mode);
			} while (!isSpinBetweenLimits(&spin[i], limits));
		}
		break;
	case GEN_FIXED_ANGLES:
		for (short i = 0; i < NUMBER_OF_BLACKHOLES; i++) {
			do {
				spin[i].magnitude = randomBetween(limits[MIN].magnitude, limits[MAX].magnitude);
				spin[i].azimuth[FIXED] = randomBetween(limits[MIN].azimuth[FIXED],
						limits[MAX].azimuth[FIXED]);
				spin[i].inclination[FIXED] = randomBetween(limits[MIN].inclination[FIXED],
						limits[MAX].inclination[FIXED]);
				convertSpin(&spin[i], inclination, mode);
			} while (!isSpinBetweenLimits(&spin[i], limits));
		}
		break;
	case GEN_PRECESSING_XYZ:
		for (short i = 0; i < NUMBER_OF_BLACKHOLES; i++) {
			do {
				for (short j = 0; j < DIMENSION; j++) {
					spin[i].component[FIXED][j] = randomBetween(limits[MIN].component[FIXED][j],
							limits[MAX].component[FIXED][j]);
				}
				convertSpin(&spin[i], inclination, mode);
			} while (!isSpinBetweenLimits(&spin[i], limits));
		}
		break;
	case GEN_PRECESSING_ANGLES:
		for (short i = 0; i < NUMBER_OF_BLACKHOLES; i++) {
			do {
				spin[i].magnitude = randomBetween(limits[MIN].magnitude, limits[MAX].magnitude);
				spin[i].azimuth[PRECESSING] = randomBetween(limits[MIN].azimuth[PRECESSING],
						limits[MAX].azimuth[PRECESSING]);
				spin[i].inclination[PRECESSING] = randomBetween(limits[MIN].inclination[PRECESSING],
						limits[MAX].inclination[PRECESSING]);
				convertSpin(&spin[i], inclination, mode);
			} while (!isSpinBetweenLimits(&spin[i], limits));
		}
		break;
	default:
		break;
	}
	SAVE_FUNCTION_FOR_TESTING();
}

void generateBinarySystemParameters(binarySystem *system, binarySystem limits[],
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
		spin[MIN] = limits[MIN].spin[0];
		spin[MAX] = limits[MAX].spin[0];
		generateSpin(&system->spin[i], spin, system->inclination, genSpin);
	}
	SAVE_FUNCTION_FOR_TESTING();
}

///@}
/// @name Printing functions
///@{

/**
 * @param file
 * @param spin
 * @param format
 */
static void printSpinParameters(FILE *file, spinParameters *spin, OutputFormat *format) {
	BACKUP_DEFINITION_LINE();
	ushort number = 3;
	char formatString[number * format->widthWithSeparator];for
(	ushort i = FIXED; i < COORDINATE_CONVENTIONS; i++) {
		setFormat(formatString, number, format);
		fprintf(file, formatString, spin->component[i][X], spin->component[i][Y],
		spin->component[i][Z]);
		fprintf(file, formatString, spin->unity[i][X], spin->unity[i][Y], spin->unity[i][Z]);
		fprintf(file, formatString, spin->azimuth[i], spin->inclination[i], spin->elevation);
		setFormatEnd(formatString, 1, format);
		fprintf(file, formatString, spin->magnitude);
	}
	SAVE_FUNCTION_FOR_TESTING();
}

void printBinarySystemParameters(FILE *file, binarySystem *system, OutputFormat *format) {
	BACKUP_DEFINITION_LINE();
	printMassParameters(file, &system->mass, format);
	printSpinParameters(file, &system->spin[0], format);
	printSpinParameters(file, &system->spin[1], format);
	ushort number = 3;
	char formatString[number * format->widthWithSeparator];setFormat
	(formatString, number, format);
	fprintf(file, formatString, system->flatness[0], system->flatness[1], system->inclination);
	setFormatEnd(formatString, number, format);
	fprintf(file, formatString, system->distance, system->coalescencePhase,
			system->coalescenceTime);
	SAVE_FUNCTION_FOR_TESTING();
}

///@}

#ifdef TEST
/// @name Testing functions
///@{

static bool isOK_magnitudeOfSpins(void) {
	double sign = 1.0;
	spinParameters spin[NUMBER_OF_BLACKHOLES];
	for (ushort i = 0; i < NUMBER_OF_BLACKHOLES; i++) {
		sign *= (i == 0 ? 1.0 : -1.0);
		for (ushort j = 0; j < COORDINATE_CONVENTIONS; j++) {
			spin[i].component[j][X] = 2.0;
			spin[i].component[j][Y] = 3.0;
			spin[i].component[j][Z] = sign * sqrt(12.0);
		}
		SAVE_FUNCTION_CALLER();
		magnitudeOfSpin(&spin[i]);
		if (spin[i].magnitude != 5.0) {
			PRINT_ERROR();
			return false;
		}
	}
	PRINT_OK();
	return true;
}

static bool isOK_isSpinBetweenLimits(void) {
	double mult = 3.0;
	spinParameters spin, limits[2];
	limits[MAX].magnitude = mult * (limits[MIN].magnitude = 0.9);
	for (ushort coordinate = 0; coordinate < COORDINATE_CONVENTIONS; coordinate++) {
		limits[MAX].azimuth[coordinate] = mult
				* (limits[MIN].azimuth[coordinate] = 0.3 + (double) coordinate);
		limits[MAX].inclination[coordinate] = mult
				* (limits[MIN].inclination[coordinate] = 0.3 + (double) coordinate);
		limits[MAX].elevation[coordinate] = mult
				* (limits[MIN].elevation[coordinate] = 0.3 + (double) coordinate);
		for (ushort dim = 0; dim < DIMENSION; dim++) {
			limits[MAX].component[coordinate][dim] = mult
					* (limits[MIN].component[coordinate][dim] = 0.3 + (double) (coordinate + dim));
			limits[MAX].unity[coordinate][dim] = mult
					* (limits[MIN].unity[coordinate][dim] = 0.1 + (double) (coordinate + dim));
		}
	}
	for (ushort i = 1; i < mult; i++) {
		spin.magnitude = limits[MAX].magnitude / (double) i;
		for (ushort coordinate = 0; coordinate < COORDINATE_CONVENTIONS; coordinate++) {
			spin.azimuth[coordinate] = limits[MAX].azimuth[coordinate] / (double) i;
			spin.inclination[coordinate] = limits[MAX].inclination[coordinate] / (double) i;
			spin.elevation[coordinate] = limits[MAX].elevation[coordinate] / (double) i;
			for (ushort dim = 0; dim < DIMENSION; dim++) {
				spin.component[coordinate][dim] = limits[MAX].component[coordinate][dim]
						/ (double) i;
				spin.unity[coordinate][dim] = limits[MAX].unity[coordinate][dim] / (double) i;
			}
		}
		SAVE_FUNCTION_CALLER();
		if (!isSpinBetweenLimits(&spin, limits)) {
			PRINT_ERROR();
			return false;
		}
	}
	double multMod = mult + 1.0;
	spin.magnitude = limits[MAX].magnitude / (double) multMod;
	for (ushort coordinate = 0; coordinate < COORDINATE_CONVENTIONS; coordinate++) {
		spin.azimuth[coordinate] = limits[MAX].azimuth[coordinate] / (double) multMod;
		spin.inclination[coordinate] = limits[MAX].inclination[coordinate] / (double) multMod;
		spin.elevation[coordinate] = limits[MAX].elevation[coordinate] / (double) multMod;
		for (ushort dim = 0; dim < DIMENSION; dim++) {
			spin.component[coordinate][dim] = limits[MAX].component[coordinate][dim]
					/ (double) multMod;
			spin.unity[coordinate][dim] = limits[MAX].unity[coordinate][dim] / (double) multMod;
		}
	}
	SAVE_FUNCTION_CALLER();
	if (isSpinBetweenLimits(&spin, limits)) {
		PRINT_ERROR();
		return false;
	}
	multMod = mult - 1.0;
	spin.magnitude = limits[MAX].magnitude * multMod;
	for (ushort coordinate = 0; coordinate < COORDINATE_CONVENTIONS; coordinate++) {
		spin.azimuth[coordinate] = limits[MAX].azimuth[coordinate] * multMod;
		spin.inclination[coordinate] = limits[MAX].inclination[coordinate] * multMod;
		spin.elevation[coordinate] = limits[MAX].elevation[coordinate] * multMod;
		for (ushort dim = 0; dim < DIMENSION; dim++) {
			spin.component[coordinate][dim] = limits[MAX].component[coordinate][dim] * multMod;
			spin.unity[coordinate][dim] = limits[MAX].unity[coordinate][dim] * multMod;
		}
	}
	SAVE_FUNCTION_CALLER();
	if (isSpinBetweenLimits(&spin, limits)) {
		PRINT_ERROR();
		return false;
	}
	PRINT_OK();
	return true;
}

bool areBinarySystemFunctionsGood(void) {
	bool isOK = true;
	if (!isOK_magnitudeOfSpins()) {
		isOK = false;
	}
	if (!isOK_isSpinBetweenLimits()) {
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
#endif	// TEST
