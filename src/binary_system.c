/**
 * @file binary_system.c
 *
 * @date Jul 20, 2011
 * @author vereb
 * @brief Binary system specific
 */

#include <math.h>
#include <assert.h>
#include "util_IO.h"
#include "util_math.h"
#include "binary_system.h"

/// @name Helper functions
///@{

/** Calculates \f$\nu, m_1, m_2\f$ from \f$\eta, M\f$.
 * @param[in,out] mass	: all mass parameters.
 */
static void m1m2ToRemainingMass(massParameters *mass) {
	SAVE_FUNCTION_FOR_TESTING();assert(mass);
	assert(mass->mass[0] > 0.0 && mass->mass[1] > 0.0);
	mass->m1_m2 = mass->mass[0] / mass->mass[1];
	mass->nu =
			mass->mass[0] < mass->mass[1] ? mass->mass[0] / mass->mass[1] :
					mass->mass[1] / mass->mass[0];
}

/** Returns the magnitude of the spins.
 * @param[in,out] spin	: spin components
 */
static void magnitudeOfSpin(spinParameters *spin) {
	SAVE_FUNCTION_FOR_TESTING();
	spin->magnitude = 0.0;
	for (short j = 0; j < DIMENSION; j++) {
		spin->magnitude += square(spin->component[FIXED][j]);
	}
	spin->magnitude = sqrt(spin->magnitude);
}

/**	Returns true, if the mass parameters are between their limits.
 * @param[in] mass		: mass parameters to examine
 * @param[in] limits	: limits of the mass parameters
 * @return true or false
 */
static bool isMassBetweenLimits(massParameters *mass, massParameters limits[]) {
	SAVE_FUNCTION_FOR_TESTING();
	bool between = true;
	for (short i = 0; i < NUMBER_OF_BLACKHOLES; i++) {
		if (limits[MIN].mass[i] > mass->mass[i] || mass->mass[i] > limits[MAX].mass[i]) {
			between = false;
		}
	}
	if (limits[MIN].eta > mass->eta || mass->eta > limits[MAX].eta) {
		between = false;
	}
	if (limits[MIN].totalMass > mass->totalMass || mass->totalMass > limits[MAX].totalMass) {
		between = false;
	}
	if (limits[MIN].chirpMass > mass->chirpMass || mass->chirpMass > limits[MAX].chirpMass) {
		between = false;
	}
	if (limits[MIN].mu > mass->mu || mass->mu > limits[MAX].mu) {
		between = false;
	}
	if (limits[MIN].nu > mass->nu || mass->nu > limits[MAX].nu) {
		between = false;
	}
	if (limits[MIN].m1_m2 > mass->m1_m2 || mass->m1_m2 > limits[MAX].m1_m2) {
		between = false;
	}
	return between;
}

/**	Returns true, if the mass parameters are between their limits.
 * @param[in] spin		: spin parameters to examine
 * @param[in] limits	: limits of the spin parameters
 * @return true or false
 */
static bool isSpinBetweenLimits(spinParameters *spin, spinParameters limits[]) {
	SAVE_FUNCTION_FOR_TESTING();
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
	return true;
}

///@}
/// @name Conversion functions
///@{

/** Converts mass parameters according
 * @param[in,out] mass		: initial and calculated mass parameters
 * @param[in]	  convert	: specifies the initial parameters
 */
static void convertMasses(massParameters *mass, conversionMode convert) {
	SAVE_FUNCTION_FOR_TESTING();assert(mass);
	switch (convert) {
	case FROM_M1M2:
		assert(mass->mass[0] > 0.0 && mass->mass[1] > 0.0);
		mass->totalMass = mass->mass[0] + mass->mass[1];
		mass->mu = mass->mass[0] * mass->mass[1] / mass->totalMass;
		mass->eta = mass->mu / mass->totalMass;
		mass->chirpMass = pow(mass->eta, 3.0 / 5.0) * mass->totalMass;
		m1m2ToRemainingMass(mass);
		break;
	case FROM_ETAM:
		assert(mass->totalMass > 0.0 && mass->eta > 0.0 && mass->eta <= 0.25);
		mass->mass[0] = (1.0 + sqrt(1.0 - 4.0 * mass->eta)) * mass->totalMass / 2.0;
		mass->mass[1] = (1.0 - sqrt(1.0 - 4.0 * mass->eta)) * mass->totalMass / 2.0;
		mass->mu = mass->eta * mass->totalMass;
		mass->chirpMass = pow(mass->eta, 3.0 / 5.0) * mass->totalMass;
		m1m2ToRemainingMass(mass);
		break;
	case FROM_ETACHIRP:
		assert(mass->chirpMass > 0.0 &&mass->eta > 0.0 && mass->eta <= 0.25);
		mass->totalMass = mass->chirpMass / pow(mass->eta, 3.0 / 5.0);
		mass->mass[0] = (1.0 + sqrt(1.0 - 4.0 * mass->eta)) * mass->totalMass / 2.0;
		mass->mass[1] = (1.0 - sqrt(1.0 - 4.0 * mass->eta)) * mass->totalMass / 2.0;
		m1m2ToRemainingMass(mass);
		break;
	default:
		break;
	}
}

/**	Converts spin components to angles in the specified frame.
 * @param[in,out]	spin	: spin parameters
 * @param[in]		index	: frame
 */
static void convertSpinFromXyzToAngles(spinParameters *spin, const ushort index) {
	SAVE_FUNCTION_FOR_TESTING();
	spin->inclination[index] = acos(spin->component[index][Z] / spin->magnitude);
	spin->azimuth[index] = -acos(
			spin->component[index][X] / spin->magnitude / sin(spin->inclination[index]));
	if (spin->component[index][Y] < 0.0) {
		spin->azimuth[index] *= -1.0;
	}
	if (!isfinite(spin->azimuth[index])) {
		spin->azimuth[index] = 0.0;
	}
}

/**	Converts spin angles to components in the specified frame.
 * @param[in,out]	spin	: spin parameters
 * @param[in]		index	: frame
 */
static void convertSpinFromAnglesToXzy(spinParameters *spin, const ushort index) {
	SAVE_FUNCTION_FOR_TESTING();
	double cosAzimuth, cosInclination;
	cosAzimuth = spin->azimuth[index] == M_PI_2 ? 0.0 : cos(spin->azimuth[index]);
	cosInclination = spin->inclination[index] == M_PI_2 ? 0.0 : cos(spin->inclination[index]);
	spin->component[index][X] = spin->magnitude * sin(spin->inclination[index]) * cosAzimuth;
	spin->component[index][Y] = spin->magnitude * sin(spin->inclination[index])
			* sin(spin->azimuth[index]);
	spin->component[index][Z] = spin->magnitude * cosInclination;
}

/**	Converts spin parameters from the fixed frame to the precessing frame.
 * @param[in,out]	spin		: spin parameters
 * @param[in]		inclination	: inclination of the precessing frame
 */
static void convertSpinFromFixedFrame(spinParameters *spin, const double inclination) {
	SAVE_FUNCTION_FOR_TESTING();assert(spin);
	double theta[2] = { +0.0, inclination };
	spin->component[PRECESSING][X] = spin->component[FIXED][X] * cos(theta[0]) * cos(theta[1])
			+ spin->component[FIXED][Y] * sin(theta[0]) * cos(theta[1])
			- spin->component[FIXED][Z] * sin(theta[1]);
	spin->component[PRECESSING][Y] = -spin->component[FIXED][X] * sin(theta[0])
			+ spin->component[FIXED][Y] * cos(theta[0]);
	spin->component[PRECESSING][Z] = spin->component[FIXED][X] * cos(theta[0]) * sin(theta[1])
			+ spin->component[FIXED][Y] * sin(theta[0]) * sin(theta[1])
			+ spin->component[FIXED][Z] * cos(theta[1]);
}

/**	Converts spin parameters from the precessing frame to the fixed frame.
 * @param[in,out]	spin		: spin parameters
 * @param[in]		inclination	: inclination of the precessing frame
 */
static void convertSpinFromPrecessionFrame(spinParameters *spin, const double inclination) {
	SAVE_FUNCTION_FOR_TESTING();assert(spin);
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
}

/**	Converts spin parameters according the conversion parameter.
 * @param[in,out]	spin		: spin parameters
 * @param[in]		inclination	: inclination fo the precessing frame
 * @param[in]		convert		: conversion mode
 */
static void convertSpin(spinParameters spin[NUMBER_OF_BLACKHOLES], const double inclination,
		conversionMode convert) {
	SAVE_FUNCTION_FOR_TESTING();assert(spin);
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
}

///@}
/// @name Generation functions
///@{

/**	Generates the mass parameters according the generation mode.
 * @param[out]	mass	: generated mass parameters
 * @param[in]	limits	: limits of the mass parameters
 * @param[in]	mode	: generation mode
 */
static void generateMass(massParameters *mass, massParameters *limits, generationMode mode) {
	SAVE_FUNCTION_FOR_TESTING();assert(mass);
	assert(limits);
	switch (mode) {
	case GEN_ETAM:
		do {
			mass->totalMass = randomBetween(limits[MIN].totalMass, limits[MAX].totalMass);
			mass->eta = randomBetween(limits[MIN].eta, limits[MAX].eta);
			convertMasses(mass, FROM_ETAM);
		} while (!isMassBetweenLimits(mass, limits));
		break;
	case GEN_M1M2:
		do {
			for (short i = 0; i < NUMBER_OF_BLACKHOLES; i++) {
				mass->mass[i] = randomBetween(limits[MIN].mass[i], limits[MAX].mass[i]);
			}
			convertMasses(mass, FROM_M1M2);
		} while (!isMassBetweenLimits(mass, limits));
		break;
	case GEN_ETACHIRP:
		do {
			mass->chirpMass = randomBetween(limits[MIN].chirpMass, limits[MAX].chirpMass);
			mass->eta = randomBetween(limits[MIN].eta, limits[MAX].eta);
			convertMasses(mass, FROM_ETACHIRP);
		} while (!isMassBetweenLimits(mass, limits));
		break;
	default:
		break;
	}
}

/** Generates the spin parameters according the generation mode.
 * @param[out]	spin		: generated spin parameters
 * @param[in]	limits		: limits of the spin parameters
 * @param[in]	inclination	: inclination of the precessing frame
 * @param[in]	mode		: generation mode
 */
static void generateSpin(spinParameters *spin, spinParameters limits[], double inclination,
		conversionMode mode) {
	SAVE_FUNCTION_FOR_TESTING();assert(spin);
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
}

void generateBinarySystemParameters(binarySystem *system, binarySystem limits[],
		generationMode genMass, generationMode genSpin) {
	SAVE_FUNCTION_FOR_TESTING();assert(system);
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
}

///@}
/// @name Printing functions
///@{

/**
 * @param file
 * @param mass
 * @param format
 */
static void printMassParameters(FILE *file, massParameters *mass, OutputFormat *format) {
	SAVE_FUNCTION_FOR_TESTING();
	ushort number = 4;
	char formatString[number * format->widthWithSeparator];
	setFormat(formatString, number, format);
	fprintf(file, formatString, mass->mass[0], mass->mass[1], mass->eta, mass->totalMass);
	fprintf(file, formatString, mass->chirpMass, mass->mu, mass->nu, mass->m1_m2);
}

/**
 * @param file
 * @param spin
 * @param format
 */
static void printSpinParameters(FILE *file, spinParameters *spin, OutputFormat *format) {
	SAVE_FUNCTION_FOR_TESTING();
	ushort number = 3;
	char formatString[number * format->widthWithSeparator];setFormat
	(formatString, number, format);
	for (ushort i = FIXED; i < COORDINATE_CONVENTIONS; i++) {
		fprintf(file, formatString, spin->component[i][X], spin->component[i][Y],
				spin->component[i][Z]);
		fprintf(file, formatString, spin->unity[i][X], spin->unity[i][Y], spin->unity[i][Z]);
		fprintf(file, formatString, spin->azimuth[i], spin->inclination[i], spin->elevation);
		fprintf(file, format->oneNumber, spin->magnitude);
	}
}

void printBinarySystemParameters(FILE *file, binarySystem *system, OutputFormat *format) {
	SAVE_FUNCTION_FOR_TESTING();
	printMassParameters(file, &system->mass, format);
	fputs("\n", file);
	printSpinParameters(file, &system->spin[0], format);
	fputs("\n", file);
	printSpinParameters(file, &system->spin[1], format);
	fputs("\n", file);
	ushort number = 3;
	char formatString[number * format->widthWithSeparator];setFormat
	(formatString, number, format);
	fprintf(file, formatString, system->flatness[0], system->flatness[1], system->inclination);
	fprintf(file, formatString, system->distance, system->coalescencePhase,
			system->coalescenceTime);
	fputs("\n", file);
}

///@}
/// @name Testing functions
///@{

static bool isOK_m1m2ToRemainingMass(void) {
	massParameters mass;
	mass.mass[0] = 1.0;
	mass.mass[1] = 2.0;
	SAVE_FUNCTION_CALLER();
	m1m2ToRemainingMass(&mass);
	if (mass.m1_m2 != 0.5) {
		PRINT_ERROR();
		return false;
	}
	if (mass.nu != 0.5) {
		PRINT_ERROR();
		return false;
	}
	mass.mass[0] = 2.0;
	mass.mass[1] = 1.0;
	SAVE_FUNCTION_CALLER();
	m1m2ToRemainingMass(&mass);
	if (mass.m1_m2 != 2.0) {
		PRINT_ERROR();
		return false;
	}
	if (mass.nu != 0.5) {
		PRINT_ERROR();
		return false;
	}
	PRINT_OK();
	return true;
}

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

static bool isOK_isMassBetweenLimits(void) {
	double mult = 3.0;
	massParameters mass, limits[2];
	limits[MAX].mass[0] = mult * (limits[MIN].mass[0] = 7.0);
	limits[MAX].mass[1] = mult * (limits[MIN].mass[1] = 3.0);
	limits[MAX].eta = mult * (limits[MIN].eta = 0.1);
	limits[MAX].totalMass = mult * (limits[MIN].totalMass = 5.0);
	limits[MAX].chirpMass = mult * (limits[MIN].chirpMass = 4.0);
	limits[MAX].mu = mult * (limits[MIN].mu = 0.4);
	limits[MAX].nu = mult * (limits[MIN].nu = 0.7);
	limits[MAX].m1_m2 = mult * (limits[MIN].m1_m2 = 1.0);
	for (ushort i = 1; i < mult; i++) {
		mass.mass[0] = limits[MAX].mass[0] / (double) i;
		mass.mass[1] = limits[MAX].mass[1] / (double) i;
		mass.eta = limits[MAX].eta / (double) i;
		mass.totalMass = limits[MAX].totalMass / (double) i;
		mass.chirpMass = limits[MAX].chirpMass / (double) i;
		mass.mu = limits[MAX].mu / (double) i;
		mass.nu = limits[MAX].nu / (double) i;
		mass.m1_m2 = limits[MAX].m1_m2 / (double) i;
		SAVE_FUNCTION_CALLER();
		if (!isMassBetweenLimits(&mass, limits)) {
			PRINT_ERROR();
			return false;
		}
	}
	double multMod = mult + 1.0;
	mass.mass[0] = limits[MAX].mass[0] / multMod;
	mass.mass[1] = limits[MAX].mass[1] / multMod;
	mass.eta = limits[MAX].eta / multMod;
	mass.totalMass = limits[MAX].totalMass / multMod;
	mass.chirpMass = limits[MAX].chirpMass / multMod;
	mass.mu = limits[MAX].mu / multMod;
	mass.nu = limits[MAX].nu / multMod;
	mass.m1_m2 = limits[MAX].m1_m2 / multMod;
	SAVE_FUNCTION_CALLER();
	if (isMassBetweenLimits(&mass, limits)) {
		PRINT_ERROR();
		return false;
	}
	multMod = mult - 1.0;
	mass.mass[0] = limits[MAX].mass[0] * multMod;
	mass.mass[1] = limits[MAX].mass[1] * multMod;
	mass.eta = limits[MAX].eta * multMod;
	mass.totalMass = limits[MAX].totalMass * multMod;
	mass.chirpMass = limits[MAX].chirpMass * multMod;
	mass.mu = limits[MAX].mu * multMod;
	mass.nu = limits[MAX].nu * multMod;
	mass.m1_m2 = limits[MAX].m1_m2 * multMod;
	SAVE_FUNCTION_CALLER();
	if (isMassBetweenLimits(&mass, limits)) {
		PRINT_ERROR();
		return false;
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

/*
static bool isOK_printMassParameters(void) {
	massParameters mass = { { 30.0, 3.0 }, 33.0, 0.23, 0.3, 1.4, 0.1, 10.0 };
	//SAVE_FUNCTION_CALLER();
	//printMassParameters(stdout, &mass, defaultFormat);
	PRINT_OK();
	return true;
}
*/

bool areBinarySystemFunctionsGood(void) {
	if (isOK_m1m2ToRemainingMass() && isOK_magnitudeOfSpins() && isOK_isMassBetweenLimits()
			&& isOK_isSpinBetweenLimits()) {
		return true;
	}
	return false;
}

///@}
