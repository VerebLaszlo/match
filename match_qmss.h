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
void find_Spin_Greater_Than1(Program_Parameters *program_Parameters, System_Parameters *parameters);

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

/**
 * Done.
 * @param prog
 * @param parameters
 * @param sig
 * @return
 */
short calc_Matches_For_ParameterPair(Program_Parameters *prog, System_Parameters *parameters,
		signalStruct *sig);

/** Calculates the running time.
 * @param program_Parameters
 * @param parameters
 */
void calc_Time(Program_Parameters *program_Parameters, System_Parameters *parameters,
		short sampling);

#endif /* MATCH_QMSS_H_ */
