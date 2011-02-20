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
} errorCodes;

/**	An enum to contains the integer type constatns.
 */
typedef enum {
	H1 = 0, H1P = 0, H1C = 1, H2 = 2, H2P = 2, H2C = 3, NUM_OF_SIGNALS = 4,
} constantValues;

/**	Structure containing the signals.
 */
typedef struct {
	double *signal[NUM_OF_SIGNALS]; ///< the signals in time domain with the following order: \f$h_{1+}\f$, \f$h_{1\times}\f$, \f$h_{2+}\f$, \f$h_{2\times}\f$.
	fftw_complex *csignal[NUM_OF_SIGNALS]; ///< the signals in frequency domain with the following order: \f$\tilde h_{1+}\f$, \f$\tilde h_{1\times}\f$, \f$\tilde h_{2+}\f$, \f$\tilde h_{2\times}\f$.
	fftw_plan plan[NUM_OF_SIGNALS]; ///< FFT plans to calculate \f$\tilde h_{ij}=F\left(h_{ij}\right)\f$
	double *psd;///<d
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

/**
 * D
 */
void product(fftw_complex out[], fftw_complex left[], fftw_complex right[], double norm[], long minfr, long maxfr);

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
 * @param[in,out] s  : as input contains the signals to be orthonormalised,
 * as output contains the orthonormalised signals, if was not specified where to strore them
 * @param[in]     minfr : the starting index
 * @param[in]     maxfr : the ending index
 */
void orthonormalise(signalStruct *out, signalStruct *s, long minfr, long maxfr);

/** Calculates the typical overlap.
 * @param[in] s   : structure containing the signals
 * @param[in] minfr : starting index
 * @param[in] maxfr : ending index
 * @return the overlap
 */
double match_Typical(signalStruct *s, long minfr, long maxfr);

/**	Calculates the M_best overlap
 * @param[in] s   : the structure containing the signals
 * @param[in] prod
 * @return the best match
 */
double match_Best(signalStruct *s, double prod[NUM_OF_SIGNALS][NUM_OF_SIGNALS]);

/**	Calculates the M_minfrimaxfr overlap
 * @param[in] s   : the structure containing the signals
 * @param[out] prod
 * @return the minfrimaxfr match
 */
double match_Worst(signalStruct *s, double prod[NUM_OF_SIGNALS][NUM_OF_SIGNALS]);

/** Calculates the two overlaps.
 * In the signalStruct structure only the signals must to be given and the length of the signals
 * @param[out] best  : the M_best overlap
 * @param[out] worst : the M_minfrimaxfr overlap
 * @param[in]  in     : the structure containing the signals
 * @param[in]  minfr   : the starting index
 * @param[in]  maxfr   : the endign index
 */
void calc_Overlap_Time(double *best, double *worst, signalStruct *in,
		long minfr, long maxfr);

#endif
