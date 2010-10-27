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
#include "variables.h"
#include "detector.h"

#define PREC "% -14.8lg "
#define PREC_PL "% -10.4lg "

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

void print_Binary_System(binary_System *sys, program_Params *params,
		FILE *stream);

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

#endif /* GENERATOR_H */
