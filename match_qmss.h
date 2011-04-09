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
#include <lal/LALSQTPNWaveformInterface.h>
#include <lal/LALNoiseModelsInspiral.h>
//#include <lal/GeneratePPNInspiral.h>
#include <lal/SimulateCoherentGW.h>	// CoherentGW
//#include <lal/LIGOMetadataTables.h>		// SimInspiralTable
#include <lal/GenerateInspiral.h>
#include <lal/LALInspiralBank.h>
#include <lal/LALDatatypes.h>		// LALStatus)
#include <lal/LALInspiral.h>
#include <lal/LALStdlib.h>
#include <lal/RealFFT.h>

short is_First;///<a

typedef enum CONSTANTS{
	NOT_FOUND = 0, FOUND = 1,
} CONSTANTS;

typedef enum {
	AMP00 = 1, AMP05=2, AMP10 = 4,
} ampContribution;

/**
 * X
 */
typedef struct LALParameters {
	LALStatus status;///<a
	CoherentGW waveform[2];///<a
	SimInspiralTable injParams[2];///<a
	PPNParamStruc ppnParams;///<a
	RandomInspiralSignalIn randIn;///<a
	short shorter;///<a
	long min_Length;///<a
	long max_Length;///<a
} LALParameters;

typedef struct System_Parameters {
	binary_System system[2];///<a
	double max_Spin;///<a
	double spin_Step;///<a
	double freq_Sampling;///<a
	double freq_Initial;///<a
	double time_Sampling;///<a
	double match_Typ;///<aa
	double match_Best;///<a
	double match_Minimax;///<a
	short shorter;///<a
	long min_Length;///<a
	long max_Length;///<a
	double freq_Min;///<a
	double freq_Max;///<a
	double freq_Step;///<a
	double min_Match;///<a
	double critical_Match;///<a
	double delta_Length;///<a
	short amp_Code[2];
	char approx[2][FILENAME_MAX];///<a
	char phase[2][FILENAME_MAX];///<a
	char spin[2][FILENAME_MAX];///<a
} System_Parameters;

typedef struct Program_Parameters {
	char (*output_Directories)[FILENAME_MAX];///<a
	long number_Of_Runs;///<a
	short precision;///<a
	short precision_To_Plot;///<a
	short width_Of_Number;///<a
	short width_Of_Number_To_Plot;///<a
	char folder[FILENAME_MAX];///<a
} Program_Parameters;

void readExactParameters(FILE *file, System_Parameters *params);

void read_Program_Parameters(Program_Parameters *parameters, System_Parameters *params,
		char *file_Name);

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

void generate_Parameters(System_Parameters *parameters, binary_System *limits);

/**
 * Done.
 * @param prog
 * @param parameters
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
void increment_Spins(System_Parameters* parameters);

/**
 * Done.
 * @param prog
 * @param parameters
 * @return
 */
short calc_Matches_For_ParameterPair(Program_Parameters *prog, System_Parameters *parameters, signalStruct *sig);

/**
 * Done
 * @param lalparams
 * @param parameters
 */
void initLALParameters(LALParameters *lalparams, System_Parameters *parameters);

/**
 * Done
 * @param lalparams
 * @param sig
 */
void createPSD(LALParameters *lalparams, signalStruct *sig);

void set_Signal_From_A1A2(short index, long elem, signalStruct *sig, LALParameters *lal, double F[]);

void set_Signal_From_HPHC(short index, long elem, signalStruct *sig, LALParameters *lal, double F[]);

#endif /* MATCH_QMSS_H_ */
