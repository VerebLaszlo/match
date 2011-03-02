/**
 * @file main_New.c
 *
 * @date Feb 17, 2011
 * @author vereb
 */

#include "match_qmss.h"

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

void proba(Program_Parameters *program_Parameters, System_Parameters *parameters,
		binary_System *limits);

void generate_Parameters2(System_Parameters *parameters, binary_System *limits);

/**
 * Done.
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char *argv[]) {
	char program_Parameters_File_Name[FILE_NAME_LENGTH];
	char parameters_File_Name[FILE_NAME_LENGTH];
	if (argc != 3) {
		puts("\"file for program parameters\" \"file for the parameter limits\"!!!");
		exit(-1);
	}
	sprintf(program_Parameters_File_Name, argv[1]);
	sprintf(parameters_File_Name, argv[2]);
	Program_Parameters program_Parameters;
	binary_System limits_Of_Parameters[2];
	System_Parameters parameters;
	read_Program_Parameters(&program_Parameters, &parameters, program_Parameters_File_Name);
	read_Parameters(limits_Of_Parameters, parameters_File_Name);
	puts("Start!!");
	//proba(&program_Parameters, &parameters, limits_Of_Parameters);
	run_Algorithm(&program_Parameters, &parameters, limits_Of_Parameters);
	puts("Done!!!");
	return 0;
}

void proba(Program_Parameters *program_Parameters, System_Parameters *parameters,
		binary_System *limits) {
	assert(program_Parameters);
	assert(parameters);
	assert(limits);
	assert(program_Parameters->number_Of_Runs >= 0);
	char temp[FILE_NAME_LENGTH];
	srand(86);
	sprintf(temp, "%s", program_Parameters->folder);
	for (long i = 0; i < program_Parameters->number_Of_Runs; i++) {
		generate_Parameters2(parameters, limits);
//		if (i == 2) {
			calc_Matches_For_ParameterPair(program_Parameters, parameters);
//		}
		sprintf(program_Parameters->folder, "%s", temp);
	}
}

void generate_Parameters2(System_Parameters *parameters, binary_System *limits) {
	assert(parameters);
	assert(limits);
	gen_Parameters(&parameters->system[0], &limits[0], &limits[1], ETAM, KAPPA_PSI);
	gen_Parameters(&parameters->system[1], &limits[0], &limits[1], ETAM, KAPPA_PSI);
}
