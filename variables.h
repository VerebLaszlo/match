/*
 * @file variables.h
 * @author Veréb László
 * @date 2010.10.24.
 */

#ifndef VARIABLES_H_
#define VARIABLES_H_

#include <assert.h>

/**
 * X
 */
typedef struct {
	long index;///<a
	double time_Sampling;///<a
	double freq_Sampling;///<a
	double freq_Initial;///<a
	double freq_Final;///<a
	double periods[2];///<a
	double periodsD;///<a
	char order[50];///<a
	char spin_Cont[10];///<a
	char name[50];///<a
	double match_Typ;///<a
	double match_Best;///<a
	double match_Worst;///<a
} program_Params;

/**
 * X
 */
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
	double gmst;///<a
	double F[2]; ///< antenna functions: \f$F_+, F_\times\f$
} antenna_Func;

/**
 * X
 */
typedef struct {
	black_Hole bh[2]; ///< parameters of the BHs
	double M; ///< the total mass of the system
	double chirpM; ///< the chirp mass of the system
	double eta; ///< symmetric mass ratio of the system
	double mu;///<a
	double incl; ///< inlination of the system
	double dist; ///< distance of the system
	double coaPhase;///<a
	double coaTime;///<a
	antenna_Func F; ///< antenna functions
} binary_System;

/**
 *		Enumeration of the detectors
 */
typedef enum detector_Enum {
	LL, LH, VIRGO, GEO600, TAMA20, TAMA300, GLASGOW, ISAS100, MPQ, CIT
} detector;

/**
 * X
 */
typedef struct detector_table {
	enum detector_Enum id;///<a
	char* name;///<a
	double nx[3];///<a
	double ny[3];///<a
	double x[3];///<a
} detector_table;

#endif /* VARIABLES_H_ */
