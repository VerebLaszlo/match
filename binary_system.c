/**
 * @file binary_system.c
 *
 * @date Jul 20, 2011
 * @author vereb
 * @brief Binary system specific
 */

#include <math.h>
#include <assert.h>
#include "util_math.h"
#include "binary_system.h"

/** Calculates \f$\nu, m_1, m_2\f$ from \f$\eta, M\f$.
 * @param[in,out] mass	: all mass parameters.
 */
static void m1m2ToRemainingMass(massParameters *mass) {
	assert(mass);
	assert(mass->mass[0] > 0.0 && mass->mass[1] > 0.0);
	mass->m1_m2 = mass->mass[0] / mass->mass[1];
	mass->nu = mass->mass[0] < mass->mass[1] ? mass->mass[0] / mass->mass[1] : mass->mass[1]
			/ mass->mass[0];
}

/** Returns the magnitude of the spins.
 * @param[in,out] spin	: spin components
 */
static void magnitudeOfSpins(spinParameters spin[]) {
	for (short i = 0; i < NUMBER_OF_BLACKHOLES; i++) {
		spin[i].magnitude = 0.0;
		for (short j = 0; j < DIMENSION; j++) {
			spin[i].magnitude += square(spin[i].component[FIXED][j]);
		}
		spin[i].magnitude = sqrt(spin[i].magnitude);
	}
}

/**	Returns true, if the mass parameters are between their limits.
 * @param[in] mass		: mass parameters to examine
 * @param[in] limits	: limits of the mass parameters
 * @return true or false
 */
static bool isMassBetweenLimits(massParameters *mass, massParameters limits[]) {
	for (short i = 0; i < NUMBER_OF_BLACKHOLES; i++) {
		if (limits[MIN].mass[i] > mass->mass[i] || mass->mass[i] > limits[MAX].mass[i]) {
			return false;
		}
	}
	if (limits[MIN].eta > mass->eta || mass->eta > limits[MAX].eta) {
		return false;
	}
	if (limits[MIN].totalMass > mass->totalMass || mass->totalMass > limits[MAX].totalMass) {
		return false;
	}
	if (limits[MIN].chirpMass > mass->chirpMass || mass->chirpMass > limits[MAX].chirpMass) {
		return false;
	}
	if (limits[MIN].mu > mass->mu || mass->mu > limits[MAX].mu) {
		return false;
	}
	if (limits[MIN].nu > mass->nu || mass->nu || limits[MAX].nu) {
		return false;
	}
	if (limits[MIN].m1_m2 > mass->m1_m2 || mass->m1_m2 > limits[MAX].m1_m2) {

	}
	return true;
}

/**	Returns true, if the mass parameters are between their limits.
 * @param[in] spin		: spin parameters to examine
 * @param[in] limits	: limits of the spin parameters
 * @return true or false
 */
static bool isSpinBetweenLimits(spinParameters *spin, spinParameters limits[]) {
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

/// @name Conversion functions
///@{

/** Converts mass parameters according
 * @param[in,out] mass		: initial and calculated mass parameters
 * @param[in]	  convert	: specifies the initial parameters
 */
static void convertMasses(massParameters *mass, conversionMode convert) {
	assert(mass);
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
static void convertSpinsFromXyzToAngles(spinParameters spin[NUMBER_OF_BLACKHOLES],
		const ushort index) {
	for (short i = 0; i < NUMBER_OF_BLACKHOLES; i++) {
		spin[i].inclination[index] = acos(spin[i].component[index][Z] / spin[i].magnitude);
		spin[i].azimuth[index] = -acos(spin[i].component[index][X] / spin[i].magnitude / sin(
				spin[i].inclination[index]));
		if (spin[i].component[index][Y] < 0.0) {
			spin[i].azimuth[index] *= -1.0;
		}
		if (!isfinite(spin[i].azimuth[index])) {
			spin[i].azimuth[index] = 0.0;
		}
	}
}

/**	Converts spin angles to components in the specified frame.
 * @param[in,out]	spin	: spin parameters
 * @param[in]		index	: frame
 */
static void convertSpinFromAnglesToXzy(spinParameters spin[], const ushort index) {
	double cosAzimuth, cosInclination;
	for (short i = 0; i < NUMBER_OF_BLACKHOLES; i++) {
		cosAzimuth = spin[i].azimuth[index] == M_PI_2 ? 0.0 : cos(spin[i].azimuth[index]);
		cosInclination = spin[i].inclination[index] == M_PI_2 ? 0.0 : cos(
				spin[i].inclination[index]);
		spin[i].component[index][X] = spin[i].magnitude * sin(spin[i].inclination[index])
				* cosAzimuth;
		spin[i].component[index][Y] = spin[i].magnitude * sin(spin[i].inclination[index]) * sin(
				spin[i].azimuth[index]);
		spin[i].component[index][Z] = spin[i].magnitude * cosInclination;
	}
}

/**	Converts spin parameters from the fixed frame to the precessing frame.
 * @param[in,out]	spin		: spin parameters
 * @param[in]		inclination	: inclination of the precessing frame
 */
static void convertSpinsFromFixedFrame(spinParameters spin[NUMBER_OF_BLACKHOLES],
		const double inclination) {
	assert(spin);
	double theta[2] = { +0.0, inclination };
	for (short i = 0; i < NUMBER_OF_BLACKHOLES; i++) {
		spin[i].component[PRECESSING][X] = spin[i].component[FIXED][X] * cos(theta[0]) * cos(
				theta[1]) + spin[i].component[FIXED][Y] * sin(theta[0]) * cos(theta[1])
				- spin[i].component[FIXED][Z] * sin(theta[1]);
		spin[i].component[PRECESSING][Y] = -spin[i].component[FIXED][X] * sin(theta[0])
				+ spin[i].component[FIXED][Y] * cos(theta[0]);
		spin[i].component[PRECESSING][Z] = spin[i].component[FIXED][X] * cos(theta[0]) * sin(
				theta[1]) + spin[i].component[FIXED][Y] * sin(theta[0]) * sin(theta[1])
				+ spin[i].component[FIXED][Z] * cos(theta[1]);
	}
}

/**	Converts spin parameters from the precessing frame to the fixed frame.
 * @param[in,out]	spin		: spin parameters
 * @param[in]		inclination	: inclination of the precessing frame
 */
static void convertSpinsFromPrecessionFrame(spinParameters spin[NUMBER_OF_BLACKHOLES],
		const double inclination) {
	assert(spin);
	double theta[2] = { -0.0, inclination };
	for (short i = 0; i < NUMBER_OF_BLACKHOLES; i++) {
		spin[i].component[FIXED][X] = spin[i].magnitude * (+spin[i].component[PRECESSING][X] * cos(
				theta[0]) * cos(theta[1]) + spin[i].component[PRECESSING][Y] * sin(theta[0])
				- spin[i].component[PRECESSING][Z] * cos(theta[0]) * sin(theta[1]));
		spin[i].component[FIXED][Y] = spin[i].magnitude * (-spin[i].component[PRECESSING][X] * sin(
				theta[0]) * sin(theta[1]) + spin[i].component[PRECESSING][Y] * cos(theta[0])
				+ spin[i].component[PRECESSING][Z] * sin(theta[0]) * sin(theta[1]));
		spin[i].component[FIXED][Z] = spin[i].magnitude * (+spin[i].component[PRECESSING][X] * sin(
				theta[1]) + spin[i].component[PRECESSING][Y] * cos(theta[1]));
	}
}

/**	Converts spin parameters according the conversion parameter.
 * @param[in,out]	spin		: spin parameters
 * @param[in]		inclination	: inclination fo the precessing frame
 * @param[in]		convert		: conversion mode
 */
static void convertSpins(spinParameters spin[NUMBER_OF_BLACKHOLES], const double inclination,
		conversionMode convert) {
	assert(spin);
	switch (convert) {
	case FROM_FIXED_XYZ:
		magnitudeOfSpins(spin);
		convertSpinsFromFixedFrame(spin, inclination);
		convertSpinsFromXyzToAngles(spin, FIXED);
		convertSpinsFromXyzToAngles(spin, PRECESSING);
		break;
	case FROM_FIXED_ANGLES:
		convertSpinFromAnglesToXzy(spin, FIXED);
		convertSpinsFromFixedFrame(spin, inclination);
		convertSpinsFromXyzToAngles(spin, PRECESSING);
		break;
	case FROM_PRECESSION_XZY:
		magnitudeOfSpins(spin);
		convertSpinsFromXyzToAngles(spin, FIXED);
		convertSpinsFromPrecessionFrame(spin, inclination);
		convertSpinsFromXyzToAngles(spin, PRECESSING);
		break;
	case FROM_PRECESSION_ANGLES:
		convertSpinFromAnglesToXzy(spin, PRECESSING);
		convertSpinsFromPrecessionFrame(spin, inclination);
		convertSpinsFromXyzToAngles(spin, FIXED);
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
	assert(mass);
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
	assert(spin);
	assert(limits);
	switch (mode) {
	case GEN_FIXED_XYZ:
		do {
			for (short i = 0; i < NUMBER_OF_BLACKHOLES; i++) {
				for (short j = 0; j < DIMENSION; j++) {
					spin[i].component[FIXED][j] = randomBetween(
							limits[MIN + 2 * i].component[FIXED][j],
							limits[MAX].component[FIXED][j]);
				}
			}
			convertSpins(spin, inclination, mode);
		} while (!isSpinBetweenLimits(spin, limits));
		break;
	case GEN_FIXED_ANGLES:
		do {
			for (short i = 0; i < NUMBER_OF_BLACKHOLES; i++) {
				spin[i].magnitude = randomBetween(limits[MIN].magnitude, limits[MAX].magnitude);
				spin[i].azimuth[FIXED] = randomBetween(limits[MIN].azimuth[FIXED],
						limits[MAX].azimuth[FIXED]);
				spin[i].inclination[FIXED] = randomBetween(limits[MIN].inclination[FIXED],
						limits[MAX].inclination[FIXED]);
			}
			convertSpins(spin, inclination, mode);
		} while (!isSpinBetweenLimits(spin, limits));
		break;
	case GEN_PRECESSING_XYZ:
		do {
			for (short i = 0; i < NUMBER_OF_BLACKHOLES; i++) {
				for (short j = 0; j < DIMENSION; j++) {
					spin[i].component[FIXED][j] = randomBetween(limits[MIN].component[FIXED][j],
							limits[MAX].component[FIXED][j]);
				}
			}
			convertSpins(spin, inclination, mode);
		} while (!isSpinBetweenLimits(spin, limits));
		break;
	case GEN_PRECESSING_ANGLES:
		do {
			for (short i = 0; i < NUMBER_OF_BLACKHOLES; i++) {
				spin[i].magnitude = randomBetween(limits[MIN].magnitude, limits[MAX].magnitude);
				spin[i].azimuth[PRECESSING] = randomBetween(limits[MIN].azimuth[PRECESSING],
						limits[MAX].azimuth[PRECESSING]);
				spin[i].inclination[PRECESSING] = randomBetween(
						limits[MIN].inclination[PRECESSING], limits[MAX].inclination[PRECESSING]);
			}
			convertSpins(spin, inclination, mode);
		} while (!isSpinBetweenLimits(spin, limits));
		break;
	default:
		break;
	}
}

/**	Generates the binary systems parameters according the generation modes.
 * @param[out]	system	: generated system parameters
 * @param[in]	limits	: limits of the system parameters
 * @param[in]	genMass	: mass generation mode
 * @param[in]	genSpin	: spin generation mode
 */
static void generateBinarySystemParameters(binarySystem *system, binarySystem limits[],
		generationMode genMass, generationMode genSpin) {
	assert(system);
	assert(limits);
	system->inclination = randomBetween(limits[MIN].inclination, limits[MAX].inclination);
	system->distance = randomBetween(limits[MIN].distance, limits[MAX].distance);
	for (short i = 0; i < NUMBER_OF_BLACKHOLES; i++) {
		system->flatness[i] = randomBetween(limits[MIN].flatness[i], limits[MAX].flatness[i]);
	}
	massParameters mass[2] = { limits[MIN].mass, limits[MAX].mass };
	generateMass(&system->mass, mass, genMass);
	spinParameters spin[4] = { limits[MIN].spin[0], limits[MAX].spin[0], limits[MIN].spin[0],
			limits[MAX].spin[0] };
	generateSpin(system->spin, spin, system->inclination, genSpin);
}

///@}
