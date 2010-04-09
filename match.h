/**
 * @file match.h
 * @author László Veréb
 * @date 2010.04.08.
 */

#include <fftw3.h>

/**
 *		The structure contains everything corresponding to one detector.
 */
typedef detector_Tag {
	size_t length;
	double *t;			///< template vector
	double *s;			///< signal vector
	double *n;			///< noise vector
	fftw_complex *ct;	///< fft of the template vector
	fftw_complex *cs;	///< fft of the signal vector
	fftw_complex *cn;	///< fft of the noise vector
	fftw_plan pt;		///< fftw-plan of the template vector
	fftw_plan ps;		///< fftw-plan of the signal vector
	fftw_plan pn;		///< fftw-plan of the noise vector
} detector;

/**
 *		The function allocates memory for the detector.
 * @param[in]	length	: the length of the vecotrs
 * @return	pointer to the detector structure
 */
void multi_Malloc(detector *det, size_t length);

/**
 *		The function deallocates the memory of the detector structure.
 * @param[in]	det	: pointer to the detector
 */
void multi_Free(detector * det);

/**
 *		The function calculates the time-correlated match afterwards.
 * @param[in]	det	: vector of detector structures
 */
void calculator(detector * det, size_t length);
