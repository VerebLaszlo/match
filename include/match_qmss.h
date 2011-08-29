/**
 * @file match_qmss.h
 *
 * @date Feb 25, 2011
 * @author vereb
 */

#ifndef MATCH_QMSS_H_
#define MATCH_QMSS_H_

#include "match.h"

#include <limits.h>
#include <time.h>

#define MOD_SPIN_INDEX 0///<C
short is_First; ///<a

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
void find_Spin_Greater_Than1(ProgramParameters *program_Parameters, System_Parameters *parameters);

/**
 * Done
 * @param parameters
 * @param limits
 */
void generate_Same_Parameters(System_Parameters *parameters, binary_System *limits,
	gen_Mode_Masses mass);

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
short incrementing_Spins(ProgramParameters *prog, System_Parameters* parameters, short index);

/**
 * Done.
 * @param prog
 * @param parameters
 * @param sig
 * @return
 */
short calc_Matches_For_ParameterPair(ProgramParameters *prog, System_Parameters *parameters,
	signalStruct *sig);

/**
 * Done.
 * @param prog
 * @param parameters
 * @param sig
 * @return
 */
short generate_Waveforms_For_Difference(ProgramParameters *prog, System_Parameters *parameters,
	signalStruct *sig);

//void find_Waveform_Errors_At_Parameter(Program_Parameters *prog, System_Parameters *parameters,
//		long index);

/** Calculates the running time.
 * @param program_Parameters
 * @param parameters
 */
void calc_Time(ProgramParameters *program_Parameters, System_Parameters *parameters,
	short sampling);

#endif /* MATCH_QMSS_H_ */
