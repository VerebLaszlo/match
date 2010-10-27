/*
 * variables.h
 *
 *  Created on: 2010.10.24.
 *      Author: vereb
 */

#ifndef VARIABLES_H_
#define VARIABLES_H_

typedef struct {
	long index;
	double time_Sampling;
	double freq_Sampling;
	double freq_Initial;
	double freq_Final;
	double periods[2];
	double periodsD;
	char order[50];
	char spin_Cont[10];
	char name[50];
	double match_Typ;
	double match_TypT;
	double match_Best;
	double match_Worst;
	double match_BestT;
	double match_WorstT;
} program_Params;

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
	double mu;
	double incl; ///< inlination of the system
	double dist; ///< distance of the system
	double coaPhase;
	double coaTime;
	antenna_Func F; ///< antenna functions
} binary_System;

/**
 *		Enumeration of the detectors
 */
typedef enum detector_Enum {
	LL, LH, VIRGO, GEO600, TAMA20, TAMA300, GLASGOW, ISAS100, MPQ, CIT
} detector;

typedef struct detector_table {
    enum detector_Enum id;
    char* name;
    double nx[3],ny[3],x[3];
} detector_table;

#endif /* VARIABLES_H_ */
