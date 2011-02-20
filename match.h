/**
 * @file match.h
 * @author László Veréb
 * @date 2010.04.08.
 */
#ifndef MATCH_H
#define MATCH_H

#include <fftw3.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/GeneratePPNInspiral.h>
#include "generator.h"
#include "detector.h"

/**	sonstants for error reporting
 */
typedef enum {
	MATCH_SUCCES = 0,
///< the code for succes
} errorCodes;

/**	An enum to contains the integer type constatns.
 */
typedef enum {
	H1 = 0, H1P = 0, H1C = 1, H2 = 2, H2P = 2, H2C = 3, NUM_OF_SIGNALS = 4,
///< the number of the signals
} constantValues;

typedef enum {
	NONE = 0, TIME = 1,
} maxMethods;

/**	Structure containing the signals.
 */
typedef struct {
	double *signal[NUM_OF_SIGNALS]; ///< the signals in time domain with the following order: $h_{1+}$, $h_{1\times}$, $h_{2+}$, $h_{2\times}$.
	fftw_complex *csignal[NUM_OF_SIGNALS]; ///< the signals in frequency domain with the following order: \f$\tilde h_{1+}\f$, \f$\tilde h_{1\times}\f$, \f$\tilde h_{2+}\f$, $\tilde h_{2\times}$.
	fftw_plan plan[NUM_OF_SIGNALS]; ///< FFT plans to calculate \f$\tilde h_{ij}=F\left(h_{ij}\right)\f$
	double *psd; ///< the psd
	long size; ///< the allocated memory for the signals
} signalStruct;

/**	Allocates memory for the signals.
 * @param[in] s		: pointer to the structure containing the signals
 * @param[in] size	: the size of the allocated memory
 * @return : the succes of the memory allocation
 */
int create_Signal_Struct(signalStruct *s, long size);

/**	Deallocates the memory of the signals.
 * @param[in] s	: pointer to the structure containing the signals
 */
void destroy_Signal_Struct(signalStruct *s);

/**	Calculates the normalized scalar product in frequency domain.
 *		Calculates the normalized scalar product in frequency domain using the
 * values between the two given frequency bins. The used formulae is:
 * \f[
 *		\newcommand{\inProd}[2]{\langle#1|#2\rangle}
 * 		\inProd{s_1}{s_2}=4\Re\int_{f_minfr}^{f_maxfr}
 * 			\frac{\tilde s_1\left(f\right)\tilde s_2^*\left(f\right)}
 * 			{S_h\left(f\right)}df
 * \f]
 * @param[in] left  : left vector
 * @param[in] right : right vector
 * @param[in] norm  : normalizing vector
 * @param[in] minfr   : starting frequency bin
 * @param[in] maxfr   : ending frequency bin
 * @return	the scalar product
 */
double inner_Product(fftw_complex left[], fftw_complex right[], double norm[],
		long minfr, long maxfr);

void product(fftw_complex out[], fftw_complex left[], fftw_complex right[],
		double norm[], long minfr, long maxfr);

/**	Normalise the given signals.
 *		The function normalise the signals according to the formula
 * \f[
 * 		\newcommand{\inProd}[2]{\langle#1|#2\rangle}
 * 		\tilde e=\frac{\tilde s}{\sqrt{\inProd{s}{s}}}
 * \f]
 * You can give a signalStruct, to store the normalised signals in, or a NULL
 * pointer, to overwrite the original signals.
 *
 * @param[out]     out : stores the normalised signals or NULL
 * @param[in, out] in  : as input contains the signals to be normalised,
 * as ouput contains the normalised signals, if was not specified where to store them
 * @param[in]      minfr : the starting index
 * @param[in]      maxfr : the ending index
 */
void normalise(signalStruct *out, signalStruct *in, long minfr, long maxfr);

/** Orthogonisate the given signals.
 *		The function orthogonise the signals according to the formula
 * \f[
 * 		\newcommand{\inProd}[2]{\langle#1|#2\rangle}
 * 		\tilde e_\perp=\tilde e_\times-e_+
 * 			\frac{\inProd{e_+}{e_\times}}{\inProd{e_+}{e_+}}
 * \f]
 * You can give a signalStruct, to store the orthogonised signals in, or a NULL
 * pointer, to overwrite the original signals.
 *
 * @param[out]     out : stores the orthogonised signals or NULL
 * @param[in, out] in  : as input contains the signals to be orthogonised,
 * as output contains the orthogonised signals, if was not specified where to store them
 * @param[in]      minfr : the starting index
 * @param[in]      maxfr : the ending index
 */
void orthogonise(signalStruct *out, signalStruct *in, long minfr, long maxfr);

/** Orthonormalise the given signals.
 * 		The function orthonormalise the signals according to the formulae
 * given in the normalise() and orthogonise() functions.
 * You can give a signalStruct, to store the orthogonised signals in, or a NULL
 * pointer, to overwrite the original signals.
 *
 * @param[out]    out : stores the orthonormalised signals
 * @param[in,out] in  : as input contains the signals to be orthonormalised,
 * as output contains the orthonormalised signals, if was not specified where to strore them
 * @param[in]     minfr : the starting index
 * @param[in]     maxfr : the ending index
 */
void orthonormalise(signalStruct *out, signalStruct *s, long minfr, long maxfr);

/** Calculates the simplest overlap.
 * \f[
 * 		\newcommand{\inProd}[2]{\langle#1|#2\rangle}
 * 		s=s_+F_++s_\timesF_\times
 * 		O=\frac{inProd{s_1}{s_2}}{sqrt{inProd{s_1}{s_1}inProd{s_2}{s_2}}}
 * \f]
 * @param[in] in    : structure containing the signals
 * @param[in] minfr : starting index
 * @param[in] maxfr : ending index
 * @return the overlap
 */
double match_Simple(signalStruct *in, long minfr, long maxfr);

/** Calculates the typical overlap.
 * @param[in] s   : structure containing the signals
 * @param[in] minfr : starting index
 * @param[in] maxfr : ending index
 * @return the overlap
 */
double match_Typical(signalStruct *s, long minfr, long maxfr, maxMethods method);

/**	Calculates the M_best overlap
 * @param[in] s   : the structure containing the signals
 * @param[in] minfr : the starting point
 * @param[in] maxfr : the ending point
 * @return the best match
 */
double match_Best(signalStruct *s, double prod[NUM_OF_SIGNALS][NUM_OF_SIGNALS]);

/**	Calculates the M_minfrimaxfr overlap
 * @param[in] s   : the structure containing the signals
 * @param[in] minfr : the starting point
 * @param[in] maxfr : the ending point
 * @return the minfrimaxfr match
 */
double match_Worst(signalStruct *s, double prod[NUM_OF_SIGNALS][NUM_OF_SIGNALS]);

/** Calculates the two overlaps.
 * In the signalStruct structure only the signals must to be given and the length of the signals
 * @param[out] best  : the M_best overlap
 * @param[out] worst : the M_minfrimaxfr overlap
 * @param[in]  s     : the structure containing the signals
 * @param[in]  minfr   : the starting index
 * @param[in]  maxfr   : the endign index
 */
void calc_Overlap(double *best, double *worst, signalStruct *s, long minfr,
		long maxfr);
void calc_Overlap_Time(double *best, double *worst, signalStruct *in,
		long minfr, long maxfr);
/**	Calculates the constatns needed by the M_best match.
 * @param[out] A     : ++, +x
 * @param[out] B     : x+, xx
 * @param[out] C     : ????????????
 * @param[in]  s     : the orthonormalised signals
 * @param[in]  minfr : the starting index
 * @param[in]  maxfr : the ending index
 */
void calculate_Constants(double *A, double *B, double *C, signalStruct *s,
		long minfr, long maxfr);

/** Calculates the overlap maxfrimased over the phases.
 * @param[in] s   : structure containing the signals
 * @param[in] minfr : starting index
 * @param[in] maxfr : ending index
 * @return the overlap maxfrimased over the phase
 */
double phase_Max(signalStruct *s, long minfr, long maxfr);

// Old Version Starts

typedef struct tagMatchStruct {
	double *signal[2];
	fftw_complex *csignal[2];
	fftw_plan plan[2];
	int length;
} matchStruct;

void mallocMatchStruct(matchStruct* m, int length);

void freeMatchStruct(matchStruct *m);

/**
 *		The structure contains everything corresponding to one detector.
 */
typedef struct detector_Tag {
	long length;
	detector det;
	double *t; ///< template vector
	double *s; ///< signal vector
	double *n; ///< noise vector
	fftw_complex *ct; ///< fft of the template vector
	fftw_complex *cs; ///< fft of the signal vector
	fftw_complex *cn; ///< fft of the noise vector
	fftw_plan pt; ///< fftw-plan of the template vector
	fftw_plan ps; ///< fftw-plan of the signal vector
	fftw_plan pn; ///< fftw-plan of the noise vector
} detector_Struct;

typedef struct Parameters {
	int count;
	char **title;
	SimInspiralTable *injParams;
	PPNParamStruc *ppnParams;
} Parameters;

/**
 *		The function allocates memory for the detector.
 * @param[in]		length	: the length of the vecotrs
 * @param[in,out]	det		: pointer to the detector structure
 * @return	pointer to the detector structure
 */
void multi_Malloc(long length, detector_Struct *det);

/**
 *		The function deallocates the memory of the detector structure.
 * @param[in]	det	: pointer to the detector
 */
void multi_Free(detector_Struct *det);

/**
 *		The function calculates the time-correlated match afterwards.
 *	The detector's vector is dynamically allocated.
 * @param[in]	index	: how many times is the function called
 * @param[in]	det		: vector of detector structures
 * @param[in]	length	: the length of the vector containing the detectors
 */
void calculate_After(long index, detector_Struct *det[], long *length);

/**
 *		The function calculates the normalized scalar product in frequency domain.
 * @param[in]	left	: left vector
 * @param[in]	right	: right vector
 * @param[in]	norm	: normalizing vector
 * @param[in]	minfr		: starting index
 * @param[in]	maxfr		: ending index
 * @return	the scalar product
 */
double scalar_freq(fftw_complex left[], fftw_complex right[], double norm[],
		long minfr, long maxfr);

/**
 *		The function use the blackman-window on the given vector.
 * @param[in]	array	: vector containing the signal
 * @param[in]	length	: the length of the vector
 * @param[out] wn		: the windowed signal
 */
void blackman(double array[], long length, double wn[]);

/**
 *		The function calculates the psd of the signal with any window-function.
 * @param[in]	n		: vector containing the signal
 * @param[in]	length	: the length of the vector
 * @param[in]	dt		: the sampling frequency of the signal
 * @param[in]	winf()	: the window-function
 * @return	the signal's psd
 */
double * psd(double n[], long length, double dt, void(*winf)(double array[],
		long length, double wn[]));

/**
 *		The function calculates the cross-correlation in time-domain.
 * @param[in]	h1		: the first vector
 * @param[in]	h2		: the second vector
 * @param[out]	dest	: the cross-correlated vector, twice the original length
 * @param[in]	length	: the input vector's length
 */
void calc_Time_Corr(double h1[], double h2[], double dest[], long length);

// Old Version Ends

#endif
