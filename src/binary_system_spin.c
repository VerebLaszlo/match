/**
 * @file binary_system_spin.c
 *
 * @date Aug 5, 2011
 * @author vereb
 * @brief
 */

#ifdef TEST
#include <stdlib.h>
#endif

#include "test.h"
#include <math.h>
#include <string.h>
#include <assert.h>
#include "util_math.h"
#include "binary_system.h"

/** Returns the magnitude of the spins.
 * @param[in,out] spin	: spin components
 */
static void magnitudeOfSpin(spinParameters *spin) {
	BACKUP_DEFINITION_LINE();    //
	assert(spin);
	spin->magnitude = 0.0;
	for (short j = 0; j < DIMENSION; j++) {
		spin->magnitude += square(spin->component[FIXED][j]);
	}
	spin->magnitude = sqrt(spin->magnitude);
	SAVE_FUNCTION_FOR_TESTING();
}

static void convertSpinsToUnity(spinParameters *spin) {
	for (ushort frame = 0; frame < COORDINATE_CONVENTIONS; frame++) {
		for (ushort dim = X; dim < DIMENSION; dim++) {
			spin->unity[frame][dim] = spin->component[frame][dim] / spin->magnitude;
		}
	}
}

/**	Returns true, if the mass parameters are between their limits.
 * @param[in] spin		: spin parameters to examine
 * @param[in] limits	: limits of the spin parameters
 * @return true or false
 */
static bool isSpinBetweenLimits(spinParameters *spin, spinParameters limits[]) {
	BACKUP_DEFINITION_LINE();    //
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
	}SAVE_FUNCTION_FOR_TESTING();
	return true;
}

/**	Converts spin components to angles in the specified frame.
 * @param[in,out]	spin	: spin parameters
 * @param[in]		index	: frame
 */
static void convertSpinFromXyzToAngles(spinParameters *spin, const ushort convention) {
	BACKUP_DEFINITION_LINE();
	if (spin->component[convention][X] == 0.0 && spin->component[convention][Y] == 0.0) {
		spin->azimuth[convention] = NAN;
	} else {
		spin->azimuth[convention] = atan2(spin->component[convention][Y],
			spin->component[convention][X]);
	}
	spin->inclination[convention] = acos(spin->component[convention][Z] / spin->magnitude);
	spin->elevation[convention] = M_PI_2 - spin->inclination[convention];
	SAVE_FUNCTION_FOR_TESTING();
}

/**	Converts spin angles to components in the specified frame.
 * @param[in,out]	spin	: spin parameters
 * @param[in]		index	: frame
 */
static void convertSpinFromAnglesToXyz(spinParameters *spin, const ushort convention) {
	BACKUP_DEFINITION_LINE();
	double cosAzimuth, cosInclination;
	cosAzimuth = spin->azimuth[convention] == M_PI_2 ? 0.0 : cos(spin->azimuth[convention]);
	cosInclination =
		spin->inclination[convention] == M_PI_2 ? 0.0 : cos(spin->inclination[convention]);
	spin->component[convention][X] = spin->magnitude * sin(spin->inclination[convention])
		* cosAzimuth;
	spin->component[convention][Y] = spin->magnitude * sin(spin->inclination[convention])
		* sin(spin->azimuth[convention]);
	spin->component[convention][Z] = spin->magnitude * cosInclination;
	SAVE_FUNCTION_FOR_TESTING();
}

/**	Converts spin parameters from the fixed frame to the precessing frame.
 * @param[in,out]	spin		: spin parameters
 * @param[in]		inclination	: inclination of the precessing frame
 */
static void convertSpinFromFixedFrame(spinParameters *spin, const double inclination) {
	BACKUP_DEFINITION_LINE();    //
	assert(spin);
	double theta[2] = { +0.0, inclination };
	if (inclination == 0.0) {
		memcpy(spin->component[PRECESSING], spin->component[FIXED], sizeof(double) * DIMENSION);
	} else if (inclination == M_PI) {
		spin->component[PRECESSING][X] = -spin->component[FIXED][X];
		spin->component[PRECESSING][Y] = spin->component[FIXED][Y];
		spin->component[PRECESSING][Z] = -spin->component[FIXED][Z];
	} else if (inclination == M_PI_4) {
		spin->component[PRECESSING][X] = spin->component[FIXED][Z] * M_SQRT2;
		spin->component[PRECESSING][Y] = spin->component[FIXED][Y];
		spin->component[PRECESSING][Z] = spin->component[FIXED][X] * M_SQRT2;
	} else if (inclination == M_PI_2) {
		spin->component[PRECESSING][X] = -spin->component[FIXED][Z];
		spin->component[PRECESSING][Y] = spin->component[FIXED][Y];
		spin->component[PRECESSING][Z] = spin->component[FIXED][X];
	} else if (inclination == M_PI_2 + M_PI_4) {
		spin->component[PRECESSING][X] = -spin->component[FIXED][X] * M_SQRT2;
		spin->component[PRECESSING][Y] = spin->component[FIXED][Y];
		spin->component[PRECESSING][Z] = -spin->component[FIXED][Z] * M_SQRT2;
	} else {
		spin->component[PRECESSING][X] = spin->component[FIXED][X] * cos(theta[0]) * cos(theta[1])
			+ spin->component[FIXED][Y] * sin(theta[0]) * cos(theta[1])
			- spin->component[FIXED][Z] * sin(theta[1]);
		spin->component[PRECESSING][Y] = -spin->component[FIXED][X] * sin(theta[0])
			+ spin->component[FIXED][Y] * cos(theta[0]);
		spin->component[PRECESSING][Z] = spin->component[FIXED][X] * cos(theta[0]) * sin(theta[1])
			+ spin->component[FIXED][Y] * sin(theta[0]) * sin(theta[1])
			+ spin->component[FIXED][Z] * cos(theta[1]);
	}SAVE_FUNCTION_FOR_TESTING();
}

/**	Converts spin parameters from the precessing frame to the fixed frame.
 * @param[in,out]	spin		: spin parameters
 * @param[in]		inclination	: inclination of the precessing frame
 */
static void convertSpinFromPrecessingFrame(spinParameters *spin, const double inclination) {
	BACKUP_DEFINITION_LINE();    //
	assert(spin);
	double theta[2] = { -0.0, inclination };
	spin->component[FIXED][X] = +spin->component[PRECESSING][X] * cos(theta[0]) * cos(theta[1])
		+ spin->component[PRECESSING][Y] * sin(theta[0])
		- spin->component[PRECESSING][Z] * cos(theta[0]) * sin(theta[1]);
	spin->component[FIXED][Y] = -spin->component[PRECESSING][X] * sin(theta[0]) * sin(theta[1])
		+ spin->component[PRECESSING][Y] * cos(theta[0])
		+ spin->component[PRECESSING][Z] * sin(theta[0]) * sin(theta[1]);
	spin->component[FIXED][Z] = +spin->component[PRECESSING][X] * sin(theta[1])
		+ spin->component[PRECESSING][Y] * cos(theta[1]);
	SAVE_FUNCTION_FOR_TESTING();
}

/**	Converts spin parameters according the conversion parameter.
 * @param[in,out]	spin		: spin parameters
 * @param[in]		inclination	: inclination fo the precessing frame
 * @param[in]		convert		: conversion mode
 */
static void convertSpin(spinParameters spin[NUMBER_OF_BLACKHOLES], const double inclination,
	conversionMode convert) {
	BACKUP_DEFINITION_LINE();    //
	assert(spin);
	switch (convert) {
	case FROM_FIXED_XYZ:
		magnitudeOfSpin(spin);
		convertSpinFromFixedFrame(spin, inclination);
		convertSpinFromXyzToAngles(spin, FIXED);
		convertSpinFromXyzToAngles(spin, PRECESSING);
		break;
	case FROM_FIXED_ANGLES:
		convertSpinFromAnglesToXyz(spin, FIXED);
		convertSpinFromFixedFrame(spin, inclination);
		convertSpinFromXyzToAngles(spin, PRECESSING);
		break;
	case FROM_PRECESSION_XZY:
		magnitudeOfSpin(spin);
		convertSpinFromXyzToAngles(spin, FIXED);
		convertSpinFromPrecessingFrame(spin, inclination);
		convertSpinFromXyzToAngles(spin, PRECESSING);
		break;
	case FROM_PRECESSION_ANGLES:
		convertSpinFromAnglesToXyz(spin, PRECESSING);
		convertSpinFromPrecessingFrame(spin, inclination);
		convertSpinFromXyzToAngles(spin, FIXED);
		break;
	default:
		break;
	}
	convertSpinsToUnity(spin);
	SAVE_FUNCTION_FOR_TESTING();
}

void generateSpin(spinParameters *spin, spinParameters limits[], double inclination,
	conversionMode mode) {
	BACKUP_DEFINITION_LINE();    //
	assert(spin);
	assert(limits);
	switch (mode) {
	case GEN_FIXED_XYZ:
		do {
			for (short j = 0; j < DIMENSION; j++) {
				spin->component[FIXED][j] = randomBetween(limits[MIN].component[FIXED][j],
					limits[MAX].component[FIXED][j]);
			}
			convertSpin(spin, inclination, mode);
		} while (!isSpinBetweenLimits(spin, limits));
		break;
	case GEN_FIXED_ANGLES:
		do {
			spin->magnitude = randomBetween(limits[MIN].magnitude, limits[MAX].magnitude);
			spin->azimuth[FIXED] = randomBetween(limits[MIN].azimuth[FIXED],
				limits[MAX].azimuth[FIXED]);
			spin->inclination[FIXED] = randomBetween(limits[MIN].inclination[FIXED],
				limits[MAX].inclination[FIXED]);
			convertSpin(spin, inclination, mode);
		} while (!isSpinBetweenLimits(spin, limits));
		break;
	case GEN_PRECESSING_XYZ:
		do {
			for (short j = 0; j < DIMENSION; j++) {
				spin->component[FIXED][j] = randomBetween(limits[MIN].component[FIXED][j],
					limits[MAX].component[FIXED][j]);
			}
			convertSpin(spin, inclination, mode);
		} while (!isSpinBetweenLimits(spin, limits));
		break;
	case GEN_PRECESSING_ANGLES:
		do {
			spin->magnitude = randomBetween(limits[MIN].magnitude, limits[MAX].magnitude);
			spin->azimuth[PRECESSING] = randomBetween(limits[MIN].azimuth[PRECESSING],
				limits[MAX].azimuth[PRECESSING]);
			spin->inclination[PRECESSING] = randomBetween(limits[MIN].inclination[PRECESSING],
				limits[MAX].inclination[PRECESSING]);
			convertSpin(spin, inclination, mode);
		} while (!isSpinBetweenLimits(spin, limits));
		break;
	default:
		break;
	}SAVE_FUNCTION_FOR_TESTING();
}

void printSpinParameters(FILE *file, spinParameters *spin, OutputFormat *format) {
	BACKUP_DEFINITION_LINE();
	ushort number = 3;
	ushort formatLength = number * format->widthWithSeparator;
	char formatString[formatLength];
	for (ushort i = FIXED; i < COORDINATE_CONVENTIONS; i++) {
		setFormat(formatString, number, format);
		fprintf(file, formatString, spin->component[i][X], spin->component[i][Y],
			spin->component[i][Z]);
		fprintf(file, formatString, spin->unity[i][X], spin->unity[i][Y], spin->unity[i][Z]);
		fprintf(file, formatString, spin->azimuth[i], spin->inclination[i], spin->elevation[i]);
		setFormatEnd(formatString, 1, format);
		fprintf(file, formatString, spin->magnitude);
	}SAVE_FUNCTION_FOR_TESTING();
}

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

static bool isOK_convertSpinFromXyzToAngles(void) {
	if (!isOK_magnitudeOfSpins()) {
		return false;
	}
	spinParameters spin;
	const ushort AZIM = 0;
	const ushort INCL = 1;
	const ushort ELEV = 2;
	double result[19][3] = { {-135.0, 90.0, 0.0}, {180.0, 135.0, -45.0}, {180.0, 90.0, 0.0}, {
			180.0, 45.0, 45.0}, {135.0, 90.0, 0.0}, {-90.0, 135.0, -45.0}, {-90.0, 90.0, 0.0}, {
			-90.0, 45.0, 45.0}, {NAN, 180.0, -90.0}, {NAN, NAN, NAN}, {NAN, 0.0, 90.0}, {90.0,
			135.0, -45.0}, {90.0, 90.0, 0.0}, {90.0, 45.0, 45.0}, {-45.0, 90.0, 0.0}, {0.0,
			135.0, -45.0}, {0.0, 90.0, 0.0}, {0.0, 45.0, 45.0}, {45.0, 90.0, 0.0}};
	memset(&spin, 0, sizeof(spinParameters));
	ushort testing = 0;
	ushort convention = FIXED;
	for (short x = -1; x <= 1; x++) {
		spin.component[convention][X] = (double) x;
		for (short y = -1; y <= 1; y++) {
			spin.component[convention][Y] = (double) y;
			for (short z = -1; z <= 1; z++) {
				spin.component[convention][Z] = (double) z;
				if (1 == abs(x) && 1 == abs(y) && 1 == abs(z)) {
					continue;
				}
				magnitudeOfSpin(&spin);
				SAVE_FUNCTION_CALLER();
				convertSpinFromXyzToAngles(&spin, convention);
				if (isnan(result[testing][AZIM])) {
					if (!isnan(spin.azimuth[convention])) {
						PRINT_ERROR();
						return false;
					}
				} else if (!isNear(radianFromDegree(result[testing][AZIM]),
						spin.azimuth[convention], EPSILON)) {
					PRINT_ERROR();
					return false;
				} else if (!isNear(radianFromDegree(result[testing][INCL]),
						spin.inclination[convention], EPSILON)) {
					PRINT_ERROR();
					return false;
				} else if (!isNear(radianFromDegree(result[testing][ELEV]),
						spin.elevation[convention], EPSILON)) {
					PRINT_ERROR();
					return false;
				}
				testing++;
			}
		}
	}
	double inc1 = acos(1.0 / sqrt(3.0));
	double ele1 = M_PI_2 - inc1;
	double inc2 = acos(-1.0 / sqrt(3.0));
	double ele2 = M_PI_2 - inc2;
	double result2[8][3] = { {-135.0, inc2, ele2}, {-135.0, inc1, ele1}, {135.0, inc2, ele2},
		{	135.0, inc1, ele1}, {-45.0, inc2, ele2}, {-45.0, inc1, ele1}, {45.0, inc2, ele2}, {
			45.0, inc1, ele1}};
	testing = 0;
	for (short x = -1; x <= 1; x++) {
		spin.component[convention][X] = (double) x;
		for (short y = -1; y <= 1; y++) {
			spin.component[convention][Y] = (double) y;
			for (short z = -1; z <= 1; z++) {
				spin.component[convention][Z] = (double) z;
				if (1 != abs(x) || 1 != abs(y) || 1 != abs(z)) {
					continue;
				}
				magnitudeOfSpin(&spin);
				SAVE_FUNCTION_CALLER();
				convertSpinFromXyzToAngles(&spin, convention);
				if (!isNear(radianFromDegree(result2[testing][AZIM]), spin.azimuth[convention],
						EPSILON)) {
					PRINT_ERROR();
					return false;
				} else if (!isNear(result2[testing][INCL], spin.inclination[convention], EPSILON)) {
					PRINT_ERROR();
					return false;
				} else if (!isNear(result2[testing][ELEV], spin.elevation[convention], EPSILON)) {
					PRINT_ERROR();
					return false;
				}
				testing++;
			}
		}
	}
	PRINT_OK();
	return true;
}

static bool isOK_convertSpinFromAnglesToXyz(void) {
	if (!isOK_convertSpinFromXyzToAngles()) {
		return false;
	}
	spinParameters spin;
	memset(&spin, 0, sizeof(spinParameters));
	ushort convention = FIXED;
	for (short x = -1; x <= 1; x++) {
		spin.component[convention][X] = (double) x;
		for (short y = -1; y <= 1; y++) {
			spin.component[convention][Y] = (double) y;
			for (short z = -1; z <= 1; z++) {
				spin.component[convention][Z] = (double) z;
				magnitudeOfSpin(&spin);
				convertSpinFromXyzToAngles(&spin, convention);
				memset(spin.component, 0, sizeof(double) * DIMENSION * COORDINATE_CONVENTIONS);
				SAVE_FUNCTION_CALLER();
				convertSpinFromAnglesToXyz(&spin, convention);
				if (spin.component[convention][X] != (double) x
					&& spin.component[convention][Y] != (double) y
					&& spin.component[convention][Z] != (double) z) {
					PRINT_ERROR();
					return false;
				}
			}
		}
	}
	PRINT_OK();
	return true;
}

static bool isOK_convertSpinFromFixedFrame(void) {
	spinParameters spin;
	ushort testing = 0;
	double incl = 0.0;
	do {
		for (short x = -1; x <= 1; x++) {
			spin.component[FIXED][X] = (double) x;
			for (short y = -1; y <= 1; y++) {
				spin.component[FIXED][Y] = (double) y;
				for (short z = -1; z <= 1; z++) {
					spin.component[FIXED][Z] = (double) z;
					SAVE_FUNCTION_CALLER();
					convertSpinFromFixedFrame(&spin, incl);
					if (incl == 0.0) {
						if (spin.component[FIXED][X] != spin.component[PRECESSING][X]
							&& spin.component[FIXED][Y] != spin.component[PRECESSING][Y]
							&& spin.component[FIXED][Z] != spin.component[PRECESSING][Z]) {
							PRINT_ERROR();
							return false;
						}
					} else if (incl == M_PI) {
						if (spin.component[FIXED][X] != -spin.component[PRECESSING][X]
							&& spin.component[FIXED][Y] != spin.component[PRECESSING][Y]
							&& spin.component[FIXED][Z] != -spin.component[PRECESSING][Z]) {
							PRINT_ERROR();
							return false;
						}
					} else if (incl == M_PI_2) {
						if (spin.component[FIXED][X] != spin.component[PRECESSING][Z]
							&& spin.component[FIXED][Y] != spin.component[PRECESSING][Y]
							&& spin.component[FIXED][Z] != spin.component[PRECESSING][X]) {
							PRINT_ERROR();
							return false;
						}
					} else if (incl == M_PI_4) {
						if (spin.component[FIXED][X] != spin.component[PRECESSING][X] * M_SQRT2
							&& spin.component[FIXED][Y] != spin.component[PRECESSING][Y]
							&& spin.component[FIXED][Z] != spin.component[PRECESSING][Z] * M_SQRT2) {
							PRINT_ERROR();
							return false;
						}
					} else if (incl == M_PI_2 + M_PI_4) {
						if (spin.component[FIXED][X] != -spin.component[PRECESSING][X] * M_SQRT2
							&& spin.component[FIXED][Y] != spin.component[PRECESSING][Y]
							&& spin.component[FIXED][Z] != -spin.component[PRECESSING][Z] * M_SQRT2) {
							PRINT_ERROR();
							return false;
						}
					} else {
						PRINT_ERROR();
						return false;
					}
					testing++;
				}
			}
		}
		incl += M_PI_4;
	}while (incl <= M_PI);
	PRINT_OK();
	return true;
}

static bool isOK_convertSpinFromPrecessionFrame(void) {
	if (!isOK_convertSpinFromFixedFrame()) {
		return false;
	}
	spinParameters spin;
	ushort testing = 0;
	double incl = 0.0;
	do {
		for (short x = -1; x <= 1; x++) {
			spin.component[FIXED][X] = (double) x;
			for (short y = -1; y <= 1; y++) {
				spin.component[FIXED][Y] = (double) y;
				for (short z = -1; z <= 1; z++) {
					spin.component[FIXED][Z] = (double) z;
					SAVE_FUNCTION_CALLER();
					convertSpinFromFixedFrame(&spin, incl);
					spin.component[FIXED][X] = spin.component[FIXED][Y] = spin.component[FIXED][Z] =
					0.0;
					SAVE_FUNCTION_CALLER();
					convertSpinFromPrecessingFrame(&spin, incl);
					if (spin.component[FIXED][X] != (double) x
						&& spin.component[FIXED][Y] != (double) y
						&& spin.component[FIXED][Z] != (double) z) {
						PRINT_ERROR();
						return false;
					}
					testing++;
				}
			}
		}
		incl += M_PI_4;
	}while (incl <= M_PI);
	PRINT_OK();
	return true;
}

static bool isOK_convertSpin(void) {
	if (!isOK_magnitudeOfSpins()) {
		return false;
	}
	if (!isOK_convertSpinFromFixedFrame()) {
		return false;
	}
	if (!isOK_convertSpinFromPrecessionFrame()) {
		return false;
	}
	if (!isOK_convertSpinFromXyzToAngles()) {
		return false;
	}
	if (!isOK_convertSpinFromAnglesToXyz()) {
		return false;
	}
	spinParameters spin;
	spin.component[FIXED][X] = -1.0;
	spin.component[FIXED][Y] = +1.0;
	spin.component[FIXED][Z] = -1.0;
	double inclination = 0.0;
	SAVE_FUNCTION_CALLER();
	convertSpin(&spin, inclination, FROM_FIXED_XYZ);
	if (memcmp(spin.component[FIXED], spin.component[PRECESSING], sizeof(double) * DIMENSION)) {
		PRINT_ERROR();
		return false;
	}
	PRINT_OK();
	return true;
}

static bool isOK_generateSpin(void) {
	if (!isOK_randomBetween()) {
		return false;
	}
	if (!isOK_convertSpin()) {
		return false;
	}
	if (!isOK_isSpinBetweenLimits()) {
		return false;
	}
	double inclination = 0.0;
	spinParameters spin, limits[2];
	memset(limits, 0, 2 * sizeof(spinParameters));
	limits[MAX].magnitude = 100.0 * (limits[MIN].magnitude = 1.0);
	limits[MAX].azimuth[FIXED] = 10000.0 * (limits[MIN].azimuth[FIXED] = 0.000314);
	limits[MAX].azimuth[PRECESSING] = 10000.0 * (limits[MIN].azimuth[PRECESSING] = 0.000314);
	limits[MAX].inclination[FIXED] = 10000.0 * (limits[MIN].inclination[FIXED] = 0.000314);
	limits[MAX].inclination[PRECESSING] = 10000.0 * (limits[MIN].inclination[PRECESSING] = 0.000314);
	limits[MAX].elevation[FIXED] = M_PI + (limits[MIN].elevation[FIXED] = -M_PI_2);
	limits[MAX].elevation[PRECESSING] = M_PI + (limits[MIN].elevation[PRECESSING] = -M_PI_4);
	limits[MAX].component[FIXED][X] = 10.0 * (limits[MIN].component[FIXED][X] = 1.0);
	limits[MAX].component[FIXED][Y] = 10.0 * (limits[MIN].component[FIXED][Y] = 1.0);
	limits[MAX].component[FIXED][Z] = 10.0 * (limits[MIN].component[FIXED][Z] = 1.0);
	limits[MAX].component[PRECESSING][X] = 40.0 + (limits[MIN].component[PRECESSING][X] = -20.0);
	limits[MAX].component[PRECESSING][Y] = 40.0 + (limits[MIN].component[PRECESSING][Y] = -20.0);
	limits[MAX].component[PRECESSING][Z] = 40.0 + (limits[MIN].component[PRECESSING][Z] = -20.0);
	SAVE_FUNCTION_CALLER();
	generateSpin(&spin, limits, inclination, GEN_FIXED_XYZ);
	if (!isSpinBetweenLimits(&spin, limits)) {
		generateSpin(&spin, limits, inclination, GEN_FIXED_XYZ);
		PRINT_ERROR();
		return false;
	}
	generateSpin(&spin, limits, inclination, GEN_FIXED_XYZ);
	PRINT_OK();
	return true;
}

bool areBinarySystemSpinFunctionsGood(void) {
	bool isOK = true;
	if (!isOK_isSpinBetweenLimits()) {
		isOK = false;
	} else if (!isOK_generateSpin()) {
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
