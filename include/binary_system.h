/**
 * @file binary_system.h
 *
 * @date Jul 20, 2011
 * @author vereb
 * @brief Binary system specific
 */

#ifndef BINARY_SYSTEM_H_
#define BINARY_SYSTEM_H_

typedef enum tagConversionMode {
	FROM_M1M2 = 0,
	FROM_ETAM,
	FROM_ETACHIRP,
	MASS_CONVERSIONS,
	FROM_FIXED_XYZ = 0,
	FROM_FIXED_ANGLES,
	FROM_PRECESSION_XZY,
	FROM_PRECESSION_ANGLES,
	SPIN_CONVERSIONS,
	FROM_OTHER = 100,
} conversionMode;

typedef enum tagGenerationMode {
	GEN_M1M2 = FROM_M1M2,
	GEN_ETAM = FROM_ETAM,
	GEN_ETACHIRP = FROM_ETACHIRP,
	GEN_FIXED_XYZ = FROM_FIXED_XYZ,
	GEN_FIXED_ANGLES = FROM_FIXED_ANGLES,
	GEN_PRECESSING_XYZ = FROM_PRECESSION_XZY,
	GEN_PRECESSING_ANGLES = FROM_PRECESSION_ANGLES,
	GEN_OTHER = FROM_OTHER,
} generationMode;

typedef enum tagSystemConstants {
	X = 0,
	Y,
	Z,
	DIMENSION,
	MIN = 0,
	MAX,
	FIXED = 0,
	PRECESSING,
	COORDINATE_CONVENTIONS = 2,
	NUMBER_OF_BLACKHOLES = 2,
} systemConstants;

typedef struct tagMassParameters {
	double mass[NUMBER_OF_BLACKHOLES];
	double eta;
	double mu;
	double totalMass;
	double chirpMass;
	double nu;
	double m1_m2;
} massParameters;

typedef struct tagSpinParameters {
	double component[COORDINATE_CONVENTIONS][DIMENSION];
	double magnitude;
	double unity[DIMENSION];
	double azimuth[COORDINATE_CONVENTIONS];
	double inclination[COORDINATE_CONVENTIONS];
	double elevation[COORDINATE_CONVENTIONS];
} spinParameters;

typedef struct tagBinarySystem {
	massParameters mass;
	spinParameters spin[NUMBER_OF_BLACKHOLES];
	double flatness[NUMBER_OF_BLACKHOLES];
	double inclination;
	double distance;
	double coalescencePhase;
	double coalescenceTime;
} binarySystem;

#endif /* BINARY_SYSTEM_H_ */
