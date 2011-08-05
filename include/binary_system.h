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

/**	Contains the spin parameters.
 */
typedef struct tagSpinParameters {
	double component[COORDINATE_CONVENTIONS][DIMENSION]; ///< components in the corresponding convention
	double magnitude; ///< magnitude of the spin
	double unity[COORDINATE_CONVENTIONS][DIMENSION]; ///< unity components in the corresponding convention
	double azimuth[COORDINATE_CONVENTIONS]; ///< azimuth in the corresponding convention, \f$\in[0,2\pi)\f$
	double inclination[COORDINATE_CONVENTIONS]; ///< inclination in the corresponding convention, \f$[0,\pi]\f$
	double elevation[COORDINATE_CONVENTIONS]; ///< elevation in the corresponding convention, \f$[-\pi/2,\pi/2]\f$
} spinParameters;

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
