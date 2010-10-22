/*
 * @file generator.h
 * @author László Veréb
 * @date 2010.03.26.
 */

#ifndef GENERATOR_H
#define GENERATOR_H

#include <float.h>
#include <stdio.h>
#include <string.h>
#include "util_math.h"
#include "detector.h"

#define PREC "% -14.8lg "

typedef struct {
	double m; ///< mass in \f$M_\odot\f$
	double theta; ///< inclination in the \f$x,y,z\f$ radiation and in \f$l,k,z\f$ frame #l=y
	double varphi; ///< azimuth in the \f$x,y,z\f$ radiation frame
	double phi; ///< azimuth in the \f$l,k,z\f$ frame where \f$\hat{l}=\hat{z}\times\hat{L_N}/\sin\iota\f$
	double kappa; ///< inclination in the \f$l,m,L_N\f$ frame
	double psi; ///< azimuth in the \f$l,m,L_N\f$ frame
	double ctheta; ///< \f$\cos\left(\theta\right)\f$
	double chi_Amp; ///< amplitude of the \f$\chi\f$ dimensionless spin
	double chi[3]; ///< the componenet of the \f$\chi\f$ dimensionless spin
} black_Hole;

/**
 * cos(theta) = det_z / r
 * cos(phi) = x / (r * sin(theta))
 * psi) =
 */
typedef struct {
	double dec; ///< declination
	double pol; ///< polarisation
	double alpha; ///< right ascension
	double gmst;
	double F[2]; ///< antenna functions: \f$F_+, F_\times\f$
} antenna_Func;

typedef struct {
	black_Hole bh[2]; ///< parameters of the BHs
	double M; ///< the total mass of the system
	double chirpM; ///< the chirp mass of the system
	double eta; ///< symmetric mass ratio of the system
	double incl; ///< inlination of the system
	double dist; ///< distance of the system
	double coaPhase;
	antenna_Func F; ///< antenna functions
} binary_System;

typedef enum {
	FROM_XYZ = 0, FROM_THETA_VPHI, FROM_CTHETA_VPHI, FROM_THETA_PHI, FROM_CTHETA_PHI, FROM_KAPPA_PSI,
} conversion_Mode_Spins;

typedef enum {
	FROM_M1M2 = 0, FROM_ETAM, FROM_ETACHIRP,
} conversion_Mode_Masses;

typedef enum {
	ETAM = 0, ETACHRIP, M1M2,
} gen_Mode_Masses;

typedef enum {
	XYZ = 0, THETA_VPHI, THETA_PHI, KAPPA_PSI,
} gen_Mode_Spin;

void init_Binary_System(binary_System *min, binary_System *max);

void print_Binary_System(binary_System *sys, FILE *stream);

void convert_Spins(binary_System *sys, conversion_Mode_Spins mode);

void convert_Masses(binary_System *sys, conversion_Mode_Masses mode);

void calc_Antenna_Function(binary_System *sys);

// generator functions

void check_Borders(binary_System *min, binary_System *max);

void gen_Mass(binary_System *sys, binary_System *min, binary_System *max,
		gen_Mode_Masses mode);

void gen_Chi(binary_System *sys, binary_System *min, binary_System *max, gen_Mode_Spin);

void gen_Sys(binary_System *sys, binary_System *min, binary_System *max);

void gen_Parameters(binary_System *sys, binary_System *min, binary_System *max,
		gen_Mode_Masses mass, gen_Mode_Spin spin);

/**
 *		The function calulates the given detectors response-function. Currentli supports only three.
 *	\bug	Még be kell fejezni a többi detektorra
 * @param[in]	det	: the detectors ID from the enumeration, (LL : Livingston, LH : Hanford, VIRGO: Virgo)
 * @param[in]	dec	: the declination of the source
 * @param[in]	phi	: \todo ez mi
 * @param[in]	pol	: the sources polarization angle
 * @param[out] fp	: F+
 * @param[out] fc	: Fx
 */
void calc_Response_For_Detector(detector det, binary_System *sys);

#endif /* GENERATOR_H */
