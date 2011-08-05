/**
 * @file binary_system_util.h
 *
 * @date Aug 5, 2011
 * @author vereb
 * @brief
 */

#ifndef BINARY_SYSTEM_UTIL_H_
#define BINARY_SYSTEM_UTIL_H_


/**	Conversion mode codes.
 */
typedef enum tagConversionMode {
	FROM_M1M2 = 0, ///< convert from \f$m_1, m_2\f$
	FROM_ETAM, ///< convert from \f$\eta, M\f$
	FROM_ETACHIRP, ///< convert from \f$\eta, \mathcal{M}\f$
	MASS_CONVERSIONS, ///< number of mass conversion modes
	FROM_FIXED_XYZ = 0, ///< convert from \f$x, y, z\f$ in fixed convention
	FROM_FIXED_ANGLES, ///< convert from inclination and azimuth in fixed convention
	FROM_PRECESSION_XZY, ///< convert from \f$x, y, z\f$ in precessing convention
	FROM_PRECESSION_ANGLES, ///< convert from inclination and azimuth in precessing convention
	SPIN_CONVERSIONS, ///< number of spin conversion modes
	FROM_OTHER = 100,
///< ???
} conversionMode;

/** Generation mode codes.
 */
typedef enum tagGenerationMode {
	GEN_M1M2 = FROM_M1M2, ///< generate in \f$m_1, m_2\f$
	GEN_ETAM = FROM_ETAM, ///< generate in \f$\eta, M\f$
	GEN_ETACHIRP = FROM_ETACHIRP, ///< generate in \f$\eta, \mathcal{M}\f$
	MASS_GENERATIONS = MASS_CONVERSIONS, ///< number of mass generation modes
	GEN_FIXED_XYZ = FROM_FIXED_XYZ, ///< generate in \f$x, y, z\f$ in fixed convention
	GEN_FIXED_ANGLES = FROM_FIXED_ANGLES, ///< generate in inclination and azimuth in fixed convention
	GEN_PRECESSING_XYZ = FROM_PRECESSION_XZY, ///< generate in \f$x, y, z\f$ in precessing convention
	GEN_PRECESSING_ANGLES = FROM_PRECESSION_ANGLES, ///< generate in inclination and azimuth in precessing convention
	SPIN_GENERATIONS = SPIN_CONVERSIONS, ///< number of spin generation modes
	GEN_OTHER = FROM_OTHER,
///< ???
} generationMode;

/**	Various constants for the binary system.
 */
typedef enum tagSystemConstants {
	X = 0, ///< x component
	Y, ///< y component
	Z, ///< z component
	DIMENSION, ///< number of the dimensions
	MIN = 0, ///< MIN
	MAX, ///< MAX
	FIXED = 0, ///< fixed convention
	PRECESSING, ///< precessing convention
	COORDINATE_CONVENTIONS = 2, ///< number of the spin conventions
	NUMBER_OF_BLACKHOLES = 2,
///< number of the blackholes
///< NUMBER_OF_BLACKHOLES
} systemConstants;

#endif /* BINARY_SYSTEM_UTIL_H_ */
