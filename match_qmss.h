/**
 * @file match_qmss.h
 *
 * @date Feb 25, 2011
 * @author vereb
 */

#ifndef MATCH_QMSS_H_
#define MATCH_QMSS_H_

#include "match.h"

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

/**
 * X
 */
typedef enum Constants {
	EXTRA_CHARACTERS = 6, FILE_NAME_LENGTH = 100,
} Constants;

/**
 * X
 */
typedef struct Program_Parameters {
	char (*output_Directories)[FILE_NAME_LENGTH];///<a
	long number_Of_Runs;///<a
	short precision;///<a
	short precision_To_Plot;///<a
	short width_Of_Number;///<a
	short width_Of_Number_To_Plot;///<a
	char folder[FILE_NAME_LENGTH];///<a
} Program_Parameters;

/**
 * X
 */
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
	char approx[2][FILE_NAME_LENGTH];///<a
	char phase[2][FILE_NAME_LENGTH];///<a
	char amp[2][FILE_NAME_LENGTH];///<a
	char spin[2][FILE_NAME_LENGTH];///<a
} System_Parameters;

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

/**
 * Done.
 * @param parameters
 * @param params
 * @param file_Name
 */
void read_Program_Parameters(Program_Parameters *parameters, System_Parameters *params,
		char *file_Name);

/**
 * Done.
 * @param parameters
 * @param file_Name
 */
void read_Parameters(binary_System *parameters, char *file_Name);

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
void generate_Parameters(System_Parameters *parameters, binary_System *limits);

/**
 * Done.
 * @param prog
 * @param parameters
 * @return
 */
short incrementing_Spins(Program_Parameters *prog, System_Parameters* parameters);

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
short calc_Matches_For_ParameterPair(Program_Parameters *prog, System_Parameters *parameters);

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

/**
 * Done
 * @param prog
 * @param parameters
 * @param sig
 */
void write_Waves_To_Files(Program_Parameters *prog, System_Parameters *parameters,
		signalStruct *sig);

/**
 * Done.
 * @todo előbb megnyitni valahol, hogy megkezdjük írással a fájlt, és ne csak hozzáfüzzünk!!!
 * @param prog
 * @param parameters
 * @param sig
 * @param index
 */
void write_Wave_To_File(Program_Parameters *prog, System_Parameters *parameters, signalStruct *sig,
		short index);

#endif /* MATCH_QMSS_H_ */
