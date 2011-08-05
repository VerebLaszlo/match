/**
 * @file binary_system_mass.h
 *
 * @date Aug 5, 2011
 * @author vereb
 * @brief
 */

#ifndef BINARY_SYSTEM_MASS_H_
#define BINARY_SYSTEM_MASS_H_

#include <math.h>
#include <assert.h>
#include "util_IO.h"
#include "util_math.h"
#include "binary_system_util.h"

/**	Contains the mass parameters.
 */
typedef struct tagMassParameters {
	double mass[NUMBER_OF_BLACKHOLES]; ///< \f$m_1,m_2\f$ in \f$M_\odot\f$
	double totalMass; ///< \f$M=m_1+m_2\f$ in \f$M_\odot\f$
	double eta; ///< \f$\eta=m_1m_2/M^2\f$
	double mu; ///< \f$\mu=m_1m_2/M\f$ in \f$M_\odot\f$
	double chirpMass; ///< \f$\mathcal{M}=M\eta^{3/5}\f$ in \f$M_\odot\f$
	double nu; ///< \f$\min(m_1,m_2)/\max(m_1,m_2)\f$
	double m1_m2; ///< \f$m_1/m_2\f$
} massParameters;

/**	Generates the mass parameters according the generation mode.
 * @param[out]	mass	: generated mass parameters
 * @param[in]	limits	: limits of the mass parameters
 * @param[in]	mode	: generation mode
 */
void generateMass(massParameters *mass, massParameters *limits, generationMode mode);

/**
 * @param file
 * @param mass
 * @param format
 */
void printMassParameters(FILE *file, massParameters *mass, OutputFormat *format);

#ifdef TEST

bool areBinarySystemMassFunctionsGood(void);

#endif	//TEST
#endif /* BINARY_SYSTEM_MASS_H_ */
