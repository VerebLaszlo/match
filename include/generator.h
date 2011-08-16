/**
 * @file generator.h
 * @author László Veréb
 * @date 2010.03.26.
 */

#ifndef GENERATOR_H
#define GENERATOR_H

#include "binary_system.h"
#include "util_math.h"

typedef conversionMode conversion_Mode_Spins;
typedef conversionMode conversion_Mode_Masses;
typedef generationMode gen_Mode_Masses;
typedef generationMode gen_Mode_Spin;
#define FROM_KAPPA_PSI FROM_PRECESSION_ANGLES
#define FROM_THETA_VPHI FROM_FIXED_ANGLES
#define FROM_XYZ FROM_FIXED_XYZ
#define KAPPA_PSI GEN_PRECESSING_ANGLES
#define THETA_VPHI GEN_FIXED_ANGLES
#define XYZ GEN_FIXED_XYZ
#define ETAM GEN_ETAM
#define M1M2 GEN_M1M2
#define ETACHRIP GEN_ETACHIRP

#define FROM_THETA_PHI FROM_OTHER
#define THETA_PHI GEN_OTHER

/** Parameters of the black hole
 */
typedef struct {
	double m; ///< mass in \f$M_\odot\f$
	double theta; ///< inclination in the \f$x,y,z\f$ radiation and in \f$l,k,z\f$ frame
	double varphi; ///< azimuth in the \f$x,y,z\f$ radiation frame
	double phi; ///< azimuth in the \f$l,k,z\f$ frame where \f$\hat{l}=\hat{z}\times\hat{L_N}/\sin\iota\f$
	double kappa; ///< inclination in the \f$l,m,L_N\f$ frame
	double psi; ///< azimuth in the \f$l,m,L_N\f$ frame
	double ctheta; ///< \f$\cos\left(\theta\right)\f$
	double chi_Amp; ///< amplitude of the \f$\chi\f$ dimensionless spin
	double chi[3]; ///< the componenet of the \f$\chi\f$ dimensionless spin
} black_Hole;

#include "detector.h"

/** Parameters of the system
 */
typedef struct {
	black_Hole bh[2]; ///< parameters of the BHs
	double M; ///< the total mass of the system
	double chirpM; ///< the chirp mass of the system
	double eta; ///< symmetric mass ratio of the system
	double mu; ///<a
	double incl; ///< inlination of the system
	double dist; ///< distance of the system
	double coaPhase; ///<a
	double coaTime; ///<a
	DetectorParamters F; ///< antenna functions
} binary_System;

///
void convert_Spins(binary_System *sys, conversion_Mode_Spins mode);

///
void convert_Masses(binary_System *sys, conversion_Mode_Masses mode);

///
void gen_Parameters(binary_System *sys, binary_System *min, binary_System *max,
	gen_Mode_Masses mass, gen_Mode_Spin spin);

/** Reads the parameter limits of the binary from the file.
 */
void read_Binary_Parameter_Limits(FILE*file, binary_System *sys);

/** Prints the parameter limits of the binary.
 * @param file
 * @param sys
 */
void print_Binary_Parameter_Limits(FILE* file, binary_System*sys);

void print_Binary_Parameter_For_Plot(FILE* file, binary_System*sys);

#endif /* GENERATOR_H */
