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
#include "ioHandling.h"

#define SQR(a) ( (a)*(a) )
/**
 *		The structure contains everything corresponding to one detector.
 */
typedef struct detector_Tag {
	long length;
	detector det;
	double *t;			///< template vector
	double *s;			///< signal vector
	double *n;			///< noise vector
	fftw_complex *ct;	///< fft of the template vector
	fftw_complex *cs;	///< fft of the signal vector
	fftw_complex *cn;	///< fft of the noise vector
	fftw_plan pt;		///< fftw-plan of the template vector
	fftw_plan ps;		///< fftw-plan of the signal vector
	fftw_plan pn;		///< fftw-plan of the noise vector
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
 * @param[in]	min		: starting index
 * @param[in]	max		: ending index
 * @return	the scalar product
 */
double scalar_freq(fftw_complex left[], fftw_complex right[], double norm[], long min, long max);

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
double * psd(double n[], long length, double dt, void(*winf)(double array[], long length, double wn[]));

/**
 *		The function calculates the cross-correlation in time-domain.
 * @param[in]	h1		: the first vector
 * @param[in]	h2		: the second vector
 * @param[out]	dest	: the cross-correlated vector, twice the original length
 * @param[in]	length	: the input vector's length
 */
void calc_Time_Corr(double h1[], double h2[], double dest[], long length);


#endif
