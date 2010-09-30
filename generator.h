/*
 * @file generator.h
 * @author László Veréb
 * @date 2010.03.26.
 */

#ifndef GENERATOR_H_
#define GENERATOR_H_

#include "util_math.h"

typedef struct {
	double m;		///< mass in \f$M_\odot\f$
	double cth;		///< \f$\cos\left(\theta\right)\f$
	double phi;		///< \f$\phi\f$
	double chi_Amp;	///< amplitude of the \f$\chi\f$ dimensionless spin
	double chi[3];	///< the componenet of the \f$\chi\f$ dimensionless spin
} black_Hole;

/**
 * cos(theta) = det_z / r
 * cos(phi) = x / (r * sin(theta))
 * psi) =
 */
typedef struct {
	double dec;		///< declination
	double pol;		///< polarisation
	double phi;		///< ????????????????????
	double F[2];	///< antenna functions: \f$F_+, F_\times\f$
} antenna_Func;

typedef struct {
	black_Hole bh[2];	///< parameters of the BHs
	double M;		///< the total mass of the system
	double eta;		///< symmetric mass ratio of the system
	double incl;	///< inlination of the system
	double dist;	///< distance of the system
	antenna_Func F;	///< antenna functions
} binary_System;

typedef enum {
	ANGLE_TO_COMP = 0, COMP_TO_ANGLE = 1, ETAM_TO_M1M2 = 0, M1M2_TO_ETAM = 1
} conversion_Mode;

void convert_Angles_Components(binary_System *sys, conversion_Mode mode);

void convert_etaM_m1m2(binary_System *sys, conversion_Mode mode);

void calc_Antenna_Function(binary_System *sys);

// old starts

typedef double *const dpc;

/**
 *		Structure containing the parameters of BHs
 */
typedef struct {
	double m;	///< mass in Solar-mass
	double th;	///< \todo cos(theta)
	double ph;	///< \todo phi
	double sp;	///< amplitude of the spin
	double sx;	///< x component of the spin
	double sy;	///< y component of the spin
	double sz;	///< z component of the spin
} black_hole;

/**
 *		Structure containing parameters of binary system
 */
typedef struct {
	black_hole bh1;	///< parameters of the first BH
	black_hole bh2;	///< parameters of the second BH
	double M;		///< total mass
	double eta;		///< symmetric mass ratio
	double incl;	///< inclination of the system
	double dist;	///< distance of the system
	double dec;		///< \todo
	double pol;		///< polarization angle
	double phi;		///< \todo
} binary_system;

/**
 *		The function generates the BH masses in M-eta space. The extreme values are from global
 *	variables.
 * @param[out]	m1	: mass of the first BH
 * @param[out]	m2	: mass of the second BH
 * @param[out]	M	: total mass of the system
 * @param[out]	eta	: symmetric mass ration
 */
void gen_M_eta(dpc m1, dpc m2, dpc M, dpc eta);

/**
 *		The function generates the BH masses in m1-m2 space. The extreme values are from global
 *	variables.
 * @param[out]	m1	: mass of the first BH
 * @param[out]	m2	: mass of the second BH
 * @param[out]	M	: total mass of the system
 * @param[out]	eta	: symmetric mass ration
 */
void gen_m1_m2(dpc m1, dpc m2, dpc M, dpc eta);

/**
 *		The function generates spin parameters for the given BH. Use this function to create a BH
 *	structure with only the spin parameters inicialized before other generator function.
 * @param[in]	min	: the mininal values for the BH
 * @param[in]	max	: the maximal values for the BH
 * @return	BH structure with inicialized spin parameters
 */
black_hole gen_Spin(const black_hole * const min, const black_hole * const max);

/**
 *		The function generates the system's inclination.
 * @return	inclination
 */
double gen_Inclination(void);

/**
 *		The function generates the system's distance.
 * @return	distance
 */
double gen_Distance(void);

/**
 *		The function generates the system's declination.
 * @return	declination
 */
double gen_Declination(void);

/**
 *		The function generates the system's polarization angle.
 * @return	polarization angle
 */
double gen_Polarization(void);

/**
 *		The function generates the system's \todo mit
 * @return
 */
double gen_Phi(void);

/**
 *		The function generates the parameters of the binary system.
 * @return	the structure containing the parameters
 */
binary_system gen_Params(/*binary_system lower, binary_system upper*/);

/**
 *		The function copies the values from one structure to another.
 * @param[in]	source	: copiing this
 * @param[out]	dest	: copiing here
 */
void set(binary_system * const dest, const binary_system * const source);

// old ends

#endif /* GENERATOR_H_ */
