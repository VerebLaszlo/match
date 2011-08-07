/**
 * @file datatypes.h
 *
 *  Created on: Jul 16, 2011
 *      Author: vereb
 */

#ifndef DATATYPES_H_
#define DATATYPES_H_

#include <fftw3.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>

/** Constants.
 */
typedef enum CONSTANTS {
	NOT_FOUND = 0, FOUND = 1,
} CONSTANTS;

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
	double *product_Signal[NUM_OF_SIGNALS]; ///<a
	fftw_complex *csignal[NUM_OF_SIGNALS]; ///< the signals in frequency domain with the following order: \f$\tilde h_{1+}\f$, \f$\tilde h_{1\times}\f$, \f$\tilde h_{2+}\f$, \f$\tilde h_{2\times}\f$.
	fftw_plan plan[NUM_OF_SIGNALS]; ///< FFT plans to calculate \f$\tilde h_{ij}=F\left(h_{ij}\right)\f$
	double *psd; ///<d
	long length[2]; ///< lnegth of the signals
	long size; ///< the allocated memory for the signals
} signalStruct;

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

void calculate_H_From_HPHC(signalStruct *signal, double *antennaFunction);

#endif /* DATATYPES_H_ */
