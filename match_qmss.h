/**
 * @file match_qmss.h
 *
 * @date Feb 25, 2011
 * @author vereb
 */

#ifndef MATCH_QMSS_H_
#define MATCH_QMSS_H_

#include "match.h"
#include "parameters.h"

#include <limits.h>
#include <time.h>

#define MOD_SPIN_INDEX 0///<C
short is_First;///<a

/** Constants.
 */
typedef enum CONSTANTS {
	NOT_FOUND = 0, FOUND = 1,
} CONSTANTS;

/** Amplitude code
 */
typedef enum {
	AMP00 = 1, AMP05 = 2, AMP10 = 4,
} ampContribution;

/**
 * Done.
 * @param program_Parameters
 * @param parameters
 * @param limits
 */
void run_Algorithm(Program_Parameters *program_Parameters, System_Parameters *parameters,
		binary_System *limits);

/**
 * Done
 * @param parameters
 * @param limits
 */
void generate_Same_Parameters(System_Parameters *parameters, binary_System *limits);

/**
 * Done
 * @param parameters
 * @param limits
 */
void generate_Parameters(System_Parameters *parameters, binary_System *limits);

/**
 * Done.
 * @param prog
 * @param parameters
 * @param index
 * @return
 */
short incrementing_Spins(Program_Parameters *prog, System_Parameters* parameters, short index);

/** Done
 * @param system
 * @param step
 */
void increment_Spin_Of_Binary_System(binary_System *system, double step);

/** Done
 * @param parameters
 */
void increment_Spins(System_Parameters* parameters, double step);

/**
 * Done.
 * @param prog
 * @param parameters
 * @param sig
 * @return
 */
short calc_Matches_For_ParameterPair(Program_Parameters *prog, System_Parameters *parameters,
		signalStruct *sig);

/**
 * Done
 * @param lalparams
 * @param sig
 */
void createPSD(LALParameters *lalparams, signalStruct *sig);

/** Sets the signal from CoherentGW's A1A2
 * @param index
 * @param elem
 * @param sig
 * @param lal
 * @param F
 */
void setSignal_From_A1A2(short index, long elem, signalStruct *sig, LALParameters *lal, double F[]);

/** Sets the signal from CoherentGW's HPHC
 * @param index
 * @param elem
 * @param sig
 * @param lal
 * @param F
 */
void setSignal_From_HPHC(short index, long elem, signalStruct *sig, LALParameters *lal, double F[]);

#endif /* MATCH_QMSS_H_ */
