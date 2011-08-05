/**
 * @file binary_system_spin.h
 *
 * @date Aug 5, 2011
 * @author vereb
 * @brief
 */

#ifndef BINARY_SYSTEM_SPIN_H_
#define BINARY_SYSTEM_SPIN_H_

#include "util_IO.h"
#include "util_math.h"
#include "binary_system_util.h"

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

/** Generates the spin parameters according the generation mode.
 * @param[out]	spin		: generated spin parameters
 * @param[in]	limits		: limits of the spin parameters
 * @param[in]	inclination	: inclination of the precessing frame
 * @param[in]	mode		: generation mode
 */
void generateSpin(spinParameters *spin, spinParameters limits[], double inclination,
		conversionMode mode);
/**
 * @param file
 * @param spin
 * @param format
 */
void printSpinParameters(FILE *file, spinParameters *spin, OutputFormat *format);

bool areBinarySystemSpinFunctionsGood(void);

#endif /* BINARY_SYSTEM_SPIN_H_ */
