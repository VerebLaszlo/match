/**
 * @file binary_system.h
 *
 * @date Jul 20, 2011
 * @author vereb
 * @brief Binary system specific
 */

#ifndef BINARY_SYSTEM_H_
#define BINARY_SYSTEM_H_

#include "util_IO.h"
#include "util_math.h"

/**	Conversion mode codes.
 */
typedef enum massConversionModes {
	FROM_M1M2 = 0, ///< convert from \f$m_1, m_2\f$
	FROM_ETAM, ///< convert from \f$\eta, M\f$
	FROM_ETACHIRP, ///< convert from \f$\eta, \mathcal{M}\f$
	MASS_CONVERSIONS,
///< number of mass conversion modes
} massConversionMode;

/**	Conversion mode codes.
 */
typedef enum spinConversionModes {
	FROM_FIXED_XYZ = 0, ///< convert from \f$x, y, z\f$ in fixed convention
	FROM_FIXED_ANGLES, ///< convert from inclination and azimuth in fixed convention
	FROM_PRECESSION_XZY, ///< convert from \f$x, y, z\f$ in precessing convention
	FROM_PRECESSION_ANGLES, ///< convert from inclination and azimuth in precessing convention
	SPIN_CONVERSIONS,
///< number of spin conversion modes
} spinConversionMode;

/** Generation mode codes.
 */
typedef enum massGenerationMode {
	GEN_M1M2 = FROM_M1M2, ///< generate in \f$m_1, m_2\f$
	GEN_ETAM = FROM_ETAM, ///< generate in \f$\eta, M\f$
	GEN_ETACHIRP = FROM_ETACHIRP, ///< generate in \f$\eta, \mathcal{M}\f$
	MASS_GENERATIONS = MASS_CONVERSIONS,
///< number of mass generation modes
} massGenerationMode;

/** Generation mode codes.
 */
typedef enum spinGenerationMode {
	GEN_FIXED_XYZ = FROM_FIXED_XYZ, ///< generate in \f$x, y, z\f$ in fixed convention
	GEN_FIXED_ANGLES = FROM_FIXED_ANGLES, ///< generate in inclination and azimuth in fixed convention
	GEN_PRECESSING_XYZ = FROM_PRECESSION_XZY, ///< generate in \f$x, y, z\f$ in precessing convention
	GEN_PRECESSING_ANGLES = FROM_PRECESSION_ANGLES, ///< generate in inclination and azimuth in precessing convention
	SPIN_GENERATIONS = SPIN_CONVERSIONS,
///< number of spin generation modes
} spinGenerationMode;

/**	Various constants for the binary system.
 */
typedef enum tagSystemConstants {
	FIXED = 0, ///< fixed convention
	PRECESSING, ///< precessing convention
	COORDINATE_CONVENTIONS, ///< number of the spin conventions
	NUMBER_OF_BLACKHOLES = 2,
///< number of the blackholes
///< NUMBER_OF_BLACKHOLES
} systemConstants;

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
} BinarySystem;

/**	Generates the mass parameters according the generation mode.
 * @param[out]	mass	: generated mass parameters
 * @param[in]	limits	: limits of the mass parameters
 * @param[in]	mode	: generation mode
 */
void generateMass(massParameters *mass, massParameters *limits, massGenerationMode mode);

/**
 * @param file
 * @param mass
 * @param format
 */
void printMassParameters(FILE *file, massParameters *mass, OutputFormat *format);

/** Generates the spin parameters according the generation mode.
 * @param[out]	spin		: generated spin parameters
 * @param[in]	limits		: limits of the spin parameters
 * @param[in]	inclination	: inclination of the precessing frame
 * @param[in]	mode		: generation mode
 */
void generateSpin(spinParameters *spin, spinParameters limits[], double inclination,
	spinGenerationMode mode);
/**
 * @param file
 * @param spin
 * @param format
 */
void printSpinParameters(FILE *file, spinParameters *spin, OutputFormat *format);

/**	Generates the binary systems parameters according the generation modes.
 * @param[out]	system	: generated system parameters
 * @param[in]	limits	: limits of the system parameters
 * @param[in]	genMass	: mass generation mode
 * @param[in]	genSpin	: spin generation mode
 */
void generateBinarySystemParameters(BinarySystem *system, BinarySystem limits[],
	massGenerationMode genMass, spinGenerationMode genSpin);

/**
 * @param file
 * @param system
 * @param format
 */
void printBinarySystemParameters(FILE *file, BinarySystem *system, OutputFormat *format);

#ifdef TEST

bool areBinarySystemMassFunctionsGood(void);

bool areBinarySystemSpinFunctionsGood(void);

bool areBinarySystemFunctionsGood(void);

#endif // TEST
#endif /* BINARY_SYSTEM_H_ */
