/**
 * @file binary_system_mass.c
 *
 * @date Aug 5, 2011
 * @author vereb
 * @brief
 */

#include "test.h"
#include <math.h>
#include <assert.h>
#include "util_math.h"
#include "binary_system.h"

/** Calculates \f$\nu, m_1, m_2\f$ from \f$\eta, M\f$.
 * @param[in,out] mass	: all mass parameters.
 */
static void m1m2ToRemainingMass(massParameters *mass) {
	BACKUP_DEFINITION_LINE(); //
	assert(mass);
	assert(mass->mass[0] > 0.0 && mass->mass[1] > 0.0);
	mass->m1_m2 = mass->mass[0] / mass->mass[1];
	mass->mu = mass->eta * mass->totalMass;
	mass->nu =
		mass->mass[0] < mass->mass[1] ? mass->mass[0] / mass->mass[1] :
			mass->mass[1] / mass->mass[0];
	SAVE_FUNCTION_FOR_TESTING();
}

/**	Returns true, if the mass parameters are between their limits.
 * @param[in] mass		: mass parameters to examine
 * @param[in] limits	: limits of the mass parameters
 * @return true or false
 */
static bool isMassBetweenLimits(massParameters *mass, massParameters limits[]) {
	BACKUP_DEFINITION_LINE(); //
	assert(mass);
	assert(limits);
	bool between = true;
	for (short i = 0; i < NUMBER_OF_BLACKHOLES; i++) {
		if (limits[MIN].mass[i] > mass->mass[i] || mass->mass[i] > limits[MAX].mass[i]) {
			between = false;
		}
	}
	if (limits[MIN].eta > mass->eta || mass->eta > limits[MAX].eta) {
		between = false;
	} else if (limits[MIN].totalMass > mass->totalMass || mass->totalMass > limits[MAX].totalMass) {
		between = false;
	} else if (limits[MIN].chirpMass > mass->chirpMass || mass->chirpMass > limits[MAX].chirpMass) {
		between = false;
	} else if (limits[MIN].mu > mass->mu || mass->mu > limits[MAX].mu) {
		between = false;
	} else if (limits[MIN].nu > mass->nu || mass->nu > limits[MAX].nu) {
		between = false;
	} else if ((limits[MIN].m1_m2 < mass->m1_m2 || mass->m1_m2 < limits[MAX].m1_m2)
		&& (limits[MAX].m1_m2 < mass->m1_m2 || mass->m1_m2 < limits[MIN].m1_m2)) {
		between = false;
	}SAVE_FUNCTION_FOR_TESTING();
	return between;
}

static bool areMassParametersNear(massParameters *left, massParameters *right) {
	bool areEqual = true;
	if (!isNear(left->mass[0], right->mass[0], EPSILON)) {
		areEqual = false;
	} else if (!isNear(left->mass[0], right->mass[0], EPSILON)) {
		areEqual = false;
	} else if (!isNear(left->totalMass, right->totalMass, EPSILON)) {
		areEqual = false;
	} else if (!isNear(left->eta, right->eta, EPSILON)) {
		areEqual = false;
	} else if (!isNear(left->chirpMass, right->chirpMass, EPSILON)) {
		areEqual = false;
	} else if (!isNear(left->m1_m2, right->m1_m2, EPSILON)) {
		areEqual = false;
	} else if (!isNear(left->mu, right->mu, EPSILON)) {
		areEqual = false;
	} else if (!isNear(left->nu, right->nu, EPSILON)) {
		areEqual = false;
	}
	return areEqual;
}

static double calcChirpMass(double totalMass, double eta) {
	return pow(eta, 3.0 / 5.0) * totalMass;
}

/**	Converts mass parameters from \f$m_1, m_2\f$.
 * @param[in,out] mass : mass parameters
 */
static void convertMassesFromM1M2(massParameters *mass) {
	BACKUP_DEFINITION_LINE(); //
	assert(mass->mass[0] > 0.0 && mass->mass[1] > 0.0);
	mass->totalMass = mass->mass[0] + mass->mass[1];
	mass->eta = mass->mass[0] * mass->mass[1] / square(mass->totalMass);
	mass->chirpMass = calcChirpMass(mass->totalMass, mass->eta);
	m1m2ToRemainingMass(mass);
	SAVE_FUNCTION_FOR_TESTING();
}

/**	Converts mass parameters from \f$\eta, M\f$.
 * @param[in,out] mass : mass parameters
 */
static void convertMassesFromEtaM(massParameters *mass) {
	BACKUP_DEFINITION_LINE(); //
	assert(mass->totalMass > 0.0 && mass->eta > 0.0 && mass->eta <= 0.25);
	mass->mass[0] = (1.0 + sqrt(1.0 - 4.0 * mass->eta)) * mass->totalMass / 2.0;
	mass->mass[1] = (1.0 - sqrt(1.0 - 4.0 * mass->eta)) * mass->totalMass / 2.0;
	mass->chirpMass = calcChirpMass(mass->totalMass, mass->eta);
	m1m2ToRemainingMass(mass);
	SAVE_FUNCTION_FOR_TESTING();
}

/**	Converts mass parameters from \f$\eta, \mathcal{M}\f$.
 * @param[in,out] mass : mass parameters
 */
static void convertMassesFromEtaChirp(massParameters *mass) {
	BACKUP_DEFINITION_LINE(); //
	assert(mass->chirpMass > 0.0 && mass->eta > 0.0 && mass->eta <= 0.25);
	mass->totalMass = mass->chirpMass / pow(mass->eta, 3.0 / 5.0);
	mass->mass[0] = (1.0 + sqrt(1.0 - 4.0 * mass->eta)) * mass->totalMass / 2.0;
	mass->mass[1] = (1.0 - sqrt(1.0 - 4.0 * mass->eta)) * mass->totalMass / 2.0;
	m1m2ToRemainingMass(mass);
	SAVE_FUNCTION_FOR_TESTING();
}

/** Converts mass parameters according
 * @param[in,out] mass		: initial and calculated mass parameters
 * @param[in]	  convert	: specifies the initial parameters
 */
static void convertMasses(massParameters *mass, conversionMode convert) {
	BACKUP_DEFINITION_LINE(); //
	assert(mass);
	switch (convert) {
	case FROM_M1M2:
		convertMassesFromM1M2(mass);
		break;
	case FROM_ETAM:
		convertMassesFromEtaM(mass);
		break;
	case FROM_ETACHIRP:
		convertMassesFromEtaChirp(mass);
		break;
	default:
		break;
	}SAVE_FUNCTION_FOR_TESTING();
}

void generateMass(massParameters *mass, massParameters *limits, generationMode mode) {
	BACKUP_DEFINITION_LINE(); //
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
	}SAVE_FUNCTION_FOR_TESTING();
}

/**
 * @param file
 * @param mass
 * @param format
 */
void printMassParameters(FILE *file, massParameters *mass, OutputFormat *format) {
	BACKUP_DEFINITION_LINE();
	ushort number = 4;
	ushort length = number * format->widthWithSeparator;
	char formatString[length];
	setFormat(formatString, number, format);
	fprintf(file, formatString, mass->mass[0], mass->mass[1], mass->eta, mass->totalMass);
	setFormatEnd(formatString, number, format);
	fprintf(file, formatString, mass->chirpMass, mass->mu, mass->nu, mass->m1_m2);
	SAVE_FUNCTION_FOR_TESTING();
}

#ifdef TEST

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

static bool isOK_convertMassesFromM1M2(void) {
	if (!isOK_m1m2ToRemainingMass()) {
		return false;
	}
	massParameters mass, result;
	mass.mass[0] = result.mass[0] = 90.0;
	mass.mass[1] = result.mass[1] = 10.0;
	result.totalMass = 100.0;
	result.mu = 9.0;
	result.eta = 0.09;
	result.chirpMass = calcChirpMass(result.totalMass, result.eta);
	m1m2ToRemainingMass(&result);
	SAVE_FUNCTION_CALLER();
	convertMassesFromM1M2(&mass);
	if (memcmp(&mass, &result, sizeof(massParameters))) {
		PRINT_ERROR();
		return false;
	}
	PRINT_OK();
	return true;
}

static bool isOK_convertMassesEtaM(void) {
	if (!isOK_m1m2ToRemainingMass()) {
		return false;
	}
	massParameters mass, result;
	mass.totalMass = result.totalMass = 100.0;
	mass.eta = result.eta = 0.09;
	result.mass[0] = 90.0;
	result.mass[1] = 10.0;
	result.mu = 9.0;
	result.chirpMass = calcChirpMass(result.totalMass, result.eta);
	m1m2ToRemainingMass(&result);
	SAVE_FUNCTION_CALLER();
	convertMassesFromEtaM(&mass);
	if (!areMassParametersNear(&mass, &result)) {
		convertMassesFromEtaM(&mass);
		PRINT_ERROR();
		return false;
	}
	convertMassesFromEtaM(&mass);
	PRINT_OK();
	return true;
}

static bool isOK_convertMassesEtaChirp(void) {
	if (!isOK_m1m2ToRemainingMass()) {
		return false;
	}
	massParameters mass, result;
	result.totalMass = 100.0;
	mass.eta = result.eta = 0.09;
	mass.chirpMass = result.chirpMass = calcChirpMass(result.totalMass, result.eta);
	result.mass[0] = 90.0;
	result.mass[1] = 10.0;
	result.mu = 9.0;
	m1m2ToRemainingMass(&result);
	SAVE_FUNCTION_CALLER();
	convertMassesFromEtaChirp(&mass);
	if (!areMassParametersNear(&mass, &result)) {
		PRINT_ERROR();
		return false;
	}
	PRINT_OK();
	return true;
}

static bool isOK_convertMasses(void) {
	if (!isOK_convertMassesFromM1M2()) {
		PRINT_ERROR();
		return false;
	}
	if (!isOK_convertMassesEtaM()) {
		PRINT_ERROR();
		return false;
	}
	if (!isOK_convertMassesEtaChirp()) {
		PRINT_ERROR();
		return false;
	}
	massParameters mass, result;
	mass.totalMass = result.totalMass = 100.0;
	mass.eta = result.eta = 0.09;
	result.mass[0] = 90.0;
	result.mass[1] = 10.0;
	result.mu = 9.0;
	result.chirpMass = calcChirpMass(result.totalMass, result.eta);
	m1m2ToRemainingMass(&result);
	SAVE_FUNCTION_CALLER();
	convertMasses(&mass, FROM_ETAM);
	if (!areMassParametersNear(&mass, &result)) {
		PRINT_ERROR();
		return false;
	}
	PRINT_OK();
	return true;
}

static bool isOK_generateMass(void) {
	if (!isOK_randomBetween()) {
		return false;
	}
	if (!isOK_isMassBetweenLimits()) {
		return false;
	}
	if (!isOK_convertMasses()) {
		return false;
	}
	massParameters mass, limits[2];
	limits[MIN].mass[0] = 3.0 * (limits[MIN].mass[1] = 3.0);
	limits[MAX].mass[0] = 2.0 * (limits[MAX].mass[1] = 30.0);
	convertMassesFromM1M2(&limits[MIN]);
	convertMassesFromM1M2(&limits[MAX]);
	SAVE_FUNCTION_CALLER();
	generateMass(&mass, limits, GEN_M1M2);
	if (!isMassBetweenLimits(&mass, limits)) {
		generateMass(&mass, limits, GEN_M1M2);
		PRINT_ERROR();
		return false;
	}
	generateMass(&mass, limits, GEN_M1M2);
	PRINT_OK();
	return true;
}

/*
 static bool isOK_printMassParameters(void) {
 massParameters mass = { {30.0, 3.0}, 33.0, 0.23, 0.3, 1.4, 0.1, 10.0};
 SAVE_FUNCTION_CALLER();
 printMassParameters(stdout, &mass, defaultFormat);
 PRINT_OK();
 return true;
 }
 */

bool areBinarySystemMassFunctionsGood(void) {
	bool isOK = true;
	if (!isOK_m1m2ToRemainingMass()) {
		isOK = false;
	} else if (!isOK_isMassBetweenLimits()) {
		isOK = false;
	} else if (!isOK_generateMass()) {
		isOK = false;
	}
	if (isOK) {
		PRINT_OK_FILE();
	} else {
		PRINT_ERROR_FILE();
	}
	return isOK;
}

#endif	// TEST
