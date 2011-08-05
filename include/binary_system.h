/**
 * @file binary_system.h
 *
 * @date Jul 20, 2011
 * @author vereb
 * @brief Binary system specific
 */

#ifndef BINARY_SYSTEM_H_
#define BINARY_SYSTEM_H_

#include "binary_system_mass.h"
#include "binary_system_spin.h"

/**	Contains the parameters of the binary system.
 */
typedef struct tagBinarySystem {
	massParameters mass; ///< mass parameters
	spinParameters spin[NUMBER_OF_BLACKHOLES]; ///< spin parameters
	double flatness[NUMBER_OF_BLACKHOLES]; ///< flatness of the blackholes
	double inclination; ///< inclination of the total angular momentum
	double distance; ///< distance of the binary system from the detector in \f$Mpc\f$
	double coalescencePhase; ///< phase at the coalescence
	double coalescenceTime; ///< the length of the gravitational wave at the coalescence in \f$s\f$
} binarySystem;

/**	Generates the binary systems parameters according the generation modes.
 * @param[out]	system	: generated system parameters
 * @param[in]	limits	: limits of the system parameters
 * @param[in]	genMass	: mass generation mode
 * @param[in]	genSpin	: spin generation mode
 */
void generateBinarySystemParameters(binarySystem *system, binarySystem limits[],
		generationMode genMass, generationMode genSpin);

/**
 * @param file
 * @param system
 * @param format
 */
void printBinarySystemParameters(FILE *file, binarySystem *system, OutputFormat *format);

#ifdef TEST

bool areBinarySystemFunctionsGood(void);

#endif // TEST

#endif /* BINARY_SYSTEM_H_ */
