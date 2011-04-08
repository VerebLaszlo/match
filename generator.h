/**
 * @file generator.h
 * @author László Veréb
 * @date 2010.03.26.
 */

#ifndef GENERATOR_H
#define GENERATOR_H

#include <float.h>
#include "util_math.h"
#include "detector.h"

/*
 typedef struct Masses {
 double mass[2];	///< masses of the black-holes in \f$M_\odot\f$
 double total_Mass;	///< total mass of the system in \f$M_\odot\f$
 double chirp_Mass;	///< chirp mass of the system in \f$M_\odot\f$
 double eta;	///< symmetric mass ratio of the system
 double mu;	///< reduced mass in \f$M_\odot\f$
 double nu;	///< mass ratio
 double mass1_mass2;	///< \f$m_1/m_2\f$
 } Masses;

 typedef struct Vector {
 double elevation;	///< \f$[-\pi/2,\pi/2]\f$
 double inclination;	///< \f$[0,\pi]\f$
 double azimuth;	///< \f$[0,2\pi)\f$
 double magnitude;	///< \f$[0,1]\f$
 double component[3];
 double unity_Component[3];
 } Vector;

 typedef struct Binary_System {
 Masses masses;
 Vector spin[2];
 Vector angular_Momentum;
 //	coordinate system
 antenna_Func pattern;
 double distance;
 double coalescence_Time;
 double coalescence_Phase;
 } Binary_System;
 */

#define PREC "% -14.8lg "///<a
#define PREC_PL "% 10.4lg "///<a
/** Conversion modes for spin
 */
typedef enum {
	FROM_XYZ = 0, ///< FROM_XYZ
	FROM_THETA_VPHI, ///< FROM_THETA_VPHI
	FROM_CTHETA_VPHI,///< FROM_CTHETA_VPHI
	FROM_THETA_PHI, ///< FROM_THETA_PHI
	FROM_CTHETA_PHI, ///< FROM_CTHETA_PHI
	FROM_KAPPA_PSI,
///< FROM_KAPPA_PSI
} conversion_Mode_Spins;///<a

/** Conversion modes for mass
 */
typedef enum {
	FROM_M1M2 = 0, ///< from \f$m_1,m_2\f$
	FROM_ETAM, ///< from \f$\eta,M\f$
	FROM_ETACHIRP,
///< from \f$\eta,M_{chirp}\f$
} conversion_Mode_Masses;

/** Generation mode for mass
 */
typedef enum {
	ETAM = 0,///< in \f$\eta,M\f$
	ETACHRIP,///< in \f$\eta,M_{chirp}\f$
	M1M2,///< in \f$m_1,m_2\f$
} gen_Mode_Masses;///<a

/** Generation mode for spin
 */
typedef enum {
	XYZ = 0, ///< in \f$x,y,z\f$
	THETA_VPHI,///< in \f$\theta,\varphi\f$
	THETA_PHI,///< in \f$\theta,\phi\f$
	KAPPA_PSI,
///< in \f$\kappa\psi\f$
} gen_Mode_Spin;///<a

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

/** Parameters of the system
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

///
void convert_Spins(binary_System *sys, conversion_Mode_Spins mode);

///
void convert_Masses(binary_System *sys, conversion_Mode_Masses mode);

/**
 * X
 */
void check_Borders(binary_System *min, binary_System *max);

///
void gen_Mass(binary_System *sys, binary_System *min, binary_System *max, gen_Mode_Masses mode);

///
void gen_Chi(binary_System *sys, binary_System *min, binary_System *max, gen_Mode_Spin);

///
void gen_Sys(binary_System *sys, binary_System *min, binary_System *max);

///
void gen_Parameters(binary_System *sys, binary_System *min, binary_System *max,
		gen_Mode_Masses mass, gen_Mode_Spin spin);

/** Reads the parameter limits of the binary from the file.
 */
void read_Binary_Parameter_Limits(FILE*file, binary_System *sys);

#endif /* GENERATOR_H */
