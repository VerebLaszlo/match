/*
 * @file generator.h
 * @author László Veréb
 * @date 2010.03.26.
 */

#ifndef GENERATOR_H
#define GENERATOR_H

#include <float.h>
#include "util_math.h"

typedef struct {
	double m; ///< mass in \f$M_\odot\f$
	double cth; ///< \f$\cos\left(\theta\right)\f$
	double phi; ///< \f$\phi\f$
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
	double phi; ///< ????????????????????
	double F[2]; ///< antenna functions: \f$F_+, F_\times\f$
} antenna_Func;

typedef struct {
	black_Hole bh[2]; ///< parameters of the BHs
	double M; ///< the total mass of the system
	double eta; ///< symmetric mass ratio of the system
	double incl; ///< inlination of the system
	double dist; ///< distance of the system
	antenna_Func F; ///< antenna functions
} binary_System;

typedef enum {
	ANGLE_TO_COMP = 0,
	COMP_TO_ANGLE = 1,
	ETAM_TO_M1M2 = 0,
	M1M2_TO_ETAM = 1,
	ETAM = 0,
	M1M2 = 1
} conversion_Mode;

void convert_Angles_Components(binary_System *sys, conversion_Mode mode);

void convert_etaM_m1m2(binary_System *sys, conversion_Mode mode);

void calc_Antenna_Function(binary_System *sys);

// generator functions

void check_Borders(binary_System *min, binary_System *max);

void gen_Mass(binary_System *sys, binary_System *min, binary_System *max,
		conversion_Mode mode);

void gen_Chi(binary_System *sys, binary_System *min, binary_System *max);

void gen_Sys(binary_System *sys, binary_System *min, binary_System *max);

void gen_Parameters(binary_System *sys, binary_System *min, binary_System *max,
		conversion_Mode mode);

#endif /* GENERATOR_H */
