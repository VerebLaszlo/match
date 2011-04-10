/**
 * @file match.h
 * @author László Veréb
 * @date 2010.04.08.
 */

#ifndef MATCH_H
#define MATCH_H

#include <fftw3.h>
#include "generator.h"
#include "lal_wrapper.h"

/**	An enum to contains the integer type constatns.
 */
typedef enum {
	H1 = 0,
	H1P = 0,
	H1C = 1,
	H2 = 2,
	H2P = 2,
	H2C = 3,
	HPP = 0,
	HPC = 1,
	HCP = 2,
	HCC = 3,
	NUM_OF_SIGNALS = 4,
	RESPONSE1 = 4,
	RESPONSE2 = 5,
	NOS_WITH_DETECTOR_RESPONSE = 6,
} constantValues;

/**	Structure containing the signals.
 */
typedef struct {
	double *signal[NOS_WITH_DETECTOR_RESPONSE]; ///< the signals in time domain with the following order: \f$h_{1+}\f$, \f$h_{1\times}\f$, \f$h_{2+}\f$, \f$h_{2\times}\f$.
	double *product_Signal[NUM_OF_SIGNALS];///<a
	fftw_complex *csignal[NUM_OF_SIGNALS]; ///< the signals in frequency domain with the following order: \f$\tilde h_{1+}\f$, \f$\tilde h_{1\times}\f$, \f$\tilde h_{2+}\f$, \f$\tilde h_{2\times}\f$.
	fftw_plan plan[NUM_OF_SIGNALS]; ///< FFT plans to calculate \f$\tilde h_{ij}=F\left(h_{ij}\right)\f$
	double *psd;///<d
	long length[2];///< lnegth of the signals
	long size; ///< the allocated memory for the signals
} signalStruct;

/**
 * Formulae given.
 * \f{gather*}{
 * 	M_{best}=\max_{t_0}\left(\frac{A+B}{2}+\left[\left(\frac{A-B}{2}\right)^2+C^2\right]^{0.5}\right)^{0.5}\;
 * 	M_{minimax}=\max_{t_0}\left(\frac{A+B}{2}-\left[\left(\frac{A-B}{2}\right)^2+C^2\right]^{0.5}\right)^{0.5}\\
 *	A=\inProd{\tilde{e}_{1+}}{\tilde{e}_{2+}}^2+\inProd{\tilde{e}_{1+}}{\tilde{e}_{2\bot}}^2\\
 *	B=\inProd{\tilde{e}_{1\bot}}{\tilde{e}_{2+}}^2+\inProd{\tilde{e}_{1\bot}}{\tilde{e}_{2\bot}}^2\\
 *	C=\inProd{\tilde{e}_{1+}}{\tilde{e}_{2+}}+\inProd{\tilde{e}_{1\bot}}{\tilde{e}_{2+}}+\inProd{\tilde{e}_{1+}}{\tilde{e}_{2\bot}}+\inProd{\tilde{e}_{1\bot}}{\tilde{e}_{2\bot}}
 * \f}
 * @param in
 * @param min_Index
 * @param max_Index
 * @param typ
 * @param best
 * @param minimax
 */
void calc_Matches(signalStruct *in, long min_Index, long max_Index, double *typ, double *best,
		double *minimax);

/**	Allocates memory for the signals.
 * @param[in] s		: pointer to the structure containing the signals
 * @param[in] size	: the size of the allocated memory
 * @return : the succes of the memory allocation
 */
void create_Signal_Struct(signalStruct *s, long size);

/**	Allocates memory for the signals.
 * @param[in] s		: pointer to the structure containing the signals
 * @param[in] size	: the size of the allocated memory
 * @return : the succes of the memory allocation
 */
void create_Signal_Struct1(signalStruct *s, long size);

/**	Deallocates the memory of the signals.
 * @param[in] s	: pointer to the structure containing the signals
 */
void destroy_Signal_Struct(signalStruct *s);

/**	Deallocates the memory of the signals.
 * @param[in] s	: pointer to the structure containing the signals
 */
void destroy_Signal_Struct1(signalStruct *s);

/** Prints the two signal
 * @param file
 * @param sig
 * @param dt
 * @param width
 * @param precision
 */
void print_Two_Signals(FILE*file, signalStruct *sig, double dt, short width, short precision);

/** Prints the two signal and there difference.
 * @param file
 * @param sig
 * @param dt
 * @param width
 * @param precision
 */
void print_Two_Signals_And_Difference(FILE*file, signalStruct *sig, double dt, short width,
		short precision);

/** Prints the two signal with their \f$+,\times\f$ polarisations.
 * @param file
 * @param sig
 * @param dt
 * @param width
 * @param precision
 */
void print_Two_Signals_With_HPHC(FILE*file, signalStruct *sig, double dt, short width,
		short precision);

/** Sets the signal from CoherentGW's A1A2
 * @param i
 * @param sig
 * @param lal
 * @param F
 */
void setSignal_From_A1A2(short i, signalStruct *sig, LALParameters *lal, double F[]);

/** Sets the signal from CoherentGW's HPHC
 * @param i
 * @param sig
 * @param lal
 * @param F
 */
void setSignal_From_HPHC(short i, signalStruct *sig, LALParameters *lal, double F[]);

#endif
