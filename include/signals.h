/**
 * @file signals.h
 *
 * @date 2011.08.18.
 * @author László Veréb
 */

#ifndef DATATYPES_H_
#define DATATYPES_H_

#include <fftw3.h>

typedef enum {
	H1P, H1C, H2P, H2C, NUMBER_OF_SIGNALS_COMPONENTS, NUMBER_OF_SIGNALS = 2, H1 = 0, H2,
} SignalComponentCodes;

/**	Structure containing the signals.
 */
typedef struct {
	double *inTime[NUMBER_OF_SIGNALS]; ///< the two signals in time domain: \f$h_i=h_{i+}F_++h_{i\times}F_\times\f$
	double *componentsInTime[NUMBER_OF_SIGNALS_COMPONENTS]; ///< the components of the signals in time domain with the following order: \f$h_{1+}, h_{1\times}, h_{2+}, h_{2\times}\f$.
	double *product[NUMBER_OF_SIGNALS_COMPONENTS];
	fftw_complex *componentsInFrequency[NUMBER_OF_SIGNALS_COMPONENTS]; ///< the signals in frequency domain with the following order: \f$\tilde h_{1+}, \tilde h_{1\times}, \tilde h_{2+}, \tilde h_{2\times}\f$.
	fftw_plan plan[NUMBER_OF_SIGNALS_COMPONENTS]; ///< FFT plans to calculate \f$\tilde h_{ij}=F\left(h_{ij}\right)\f$
	double *powerSpectrumDensity; ///< the used power spectrum density
	long length[2]; ///< length of the signals
	long size; ///< the allocated memory for the signals
} SignalStruct;

/**	Creates the signal structure.
 * @param[out] signal : the allocated memory for the signals
 * @param[in]  size	  : the size of the allocated memories
 */
void createSignal(SignalStruct *signal, size_t size);

/**	Destroys the signal structure.
 * @param[in] signal :
 */
void destroySignal(SignalStruct *signal);

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
	HPP = 0,
	HPC = 1,
	HCP = 2,
	HCC = 3,
	NUM_OF_SIGNALS = 4,
	RESPONSE1 = 4,
	RESPONSE2 = 5,
	NOS_WITH_DETECTOR_RESPONSE = 6,
} constantValues;

/**	Allocates memory for the signals.
 * @param[in] s		: pointer to the structure containing the signals
 * @param[in] size	: the size of the allocated memory
 * @return : the succes of the memory allocation
 */
void create_Signal_Struct(SignalStruct *s, long size);

/**	Allocates memory for the signals.
 * @param[in] s		: pointer to the structure containing the signals
 * @param[in] size	: the size of the allocated memory
 * @return : the succes of the memory allocation
 */
void create_Signal_Struct1(SignalStruct *s, long size);

/**	Deallocates the memory of the signals.
 * @param[in] s	: pointer to the structure containing the signals
 */
void destroy_Signal_Struct(SignalStruct *s);

/**	Deallocates the memory of the signals.
 * @param[in] s	: pointer to the structure containing the signals
 */
void destroy_Signal_Struct1(SignalStruct *s);

/** Prints the two signal
 * @param file
 * @param sig
 * @param dt
 * @param width
 * @param precision
 */
void print_Two_Signals(FILE*file, SignalStruct *sig, double dt, short width, short precision);

/** Prints the two signal and there difference.
 * @param file
 * @param sig
 * @param dt
 * @param width
 * @param precision
 */
void print_Two_Signals_And_Difference(FILE*file, SignalStruct *sig, double dt, short width,
	short precision);

/** Prints the two signal with their \f$+,\times\f$ polarisations.
 * @param file
 * @param sig
 * @param dt
 * @param width
 * @param precision
 */
void print_Two_Signals_With_HPHC(FILE*file, SignalStruct *sig, double dt, short width,
	short precision);

void calculate_H_From_HPHC(SignalStruct *signal, double *antennaFunction);

#endif /* DATATYPES_H_ */
