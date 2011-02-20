/**
 * @file main_New.c
 *
 * @date Feb 17, 2011
 * @author vereb
 */

#include "generator.h"
#include "match.h"
#include "util.h"
#include "detector.h"
#include "match_Multi.h"

#include <stdio.h>

/*
 #include <lal/LALSQTPNWaveformInterface.h>
 #include <lal/LALNoiseModelsInspiral.h>
 #include <lal/GeneratePPNInspiral.h>
 #include <lal/SimulateCoherentGW.h>	// CoherentGW
 #include <lal/LIGOMetadataTables.h>		// SimInspiralTable
 #include <lal/GenerateInspiral.h>
 #include <lal/LALInspiralBank.h>
 #include <lal/LALDatatypes.h>		// LALStatus)
 #include <lal/LALInspiral.h>
 #include <lal/LALStdlib.h>
 #include <lal/RealFFT.h>
 */

static double max_Spin;
static double spin_Step;

typedef enum Constants {
	EXTRA_CHARACTERS, FILE_NAME_LENGTH = 30,
} Constants;

typedef struct Program_Parameters {
	char (*output_Directories)[FILE_NAME_LENGTH];
	long number_Of_Runs;
	short precision;
	short precision_To_Plot;
	short width_Of_Number;
	short width_Of_Number_To_Plot;
} Program_Parameters;

typedef binary_System System_Parameters;

/**
 * @param parameters
 */
void read_Program_Parameters(Program_Parameters *parameters);

void read_Parameters(System_Parameters *parameters);

void run_Algorithm(Program_Parameters *program_Parameters, System_Parameters *parameters);

void generate_Parameters(System_Parameters *parameters, System_Parameters *limits);

void incrementing_Spins(System_Parameters* parameters);

void increment_Spins(System_Parameters* parameters);

void calc_Matches_For_ParameterPair(System_Parameters *parameters);

/**
 * Main function.
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char *argv[]) {
	Program_Parameters program_Parameters;
	System_Parameters limits_Of_Parameters[4];
	read_Program_Parameters(&program_Parameters);
	read_Parameters(limits_Of_Parameters);
	run_Algorithm(&program_Parameters, limits_Of_Parameters);
	return 0;
}

/// @todo megírni
void read_Program_Parameters(Program_Parameters *parameters) {

}

/** @todo megírni
 * min_1, max_1, min_2, max_2
 * 0,	  1,	 2,		3
 * indexekbe kell beleolvasni a paraméter határokat
 * @param parameters
 */
void read_Parameters(System_Parameters *parameters) {

}

/// @todo megírni
void run_Algorithm(Program_Parameters *program_Parameters, System_Parameters *limits) {
	System_Parameters actual_Parameters[2];
	for (long i = 0; i < program_Parameters->number_Of_Runs; i++) {
		generate_Parameters(&actual_Parameters[0], &limits[0]);
		generate_Parameters(&actual_Parameters[1], &limits[2]);
		incrementing_Spins(actual_Parameters);
	}
}

/// @todo megírni
void generate_Parameters(System_Parameters *parameters, System_Parameters *limits) {

}

/// @todo megírni
void incrementing_Spins(System_Parameters* parameters) {
	for (; parameters->bh[0].chi_Amp < max_Spin; increment_Spins(parameters)) {
		calc_Matches_For_ParameterPair(parameters);
	}
}

void increment_Spins(System_Parameters* parameters) {
	parameters[1].bh[0].chi_Amp += spin_Step;
	parameters[1].bh[1].chi_Amp += spin_Step;
}

/// @todo megírni
void calc_Matches_For_ParameterPair(System_Parameters *parameters) {

}
