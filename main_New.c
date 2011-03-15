/**
 * @file main_New.c
 *
 * @date Feb 17, 2011
 * @author vereb
 */

#include "match_qmss.h"

#define MOD_SPIN_INDEX 0///<C
short is_First;///<a

void proba(Program_Parameters *program_Parameters, System_Parameters *parameters,
		binary_System *limits);

void generate_Parameters2(System_Parameters *parameters, binary_System *limits);

void readExactParameters(FILE *file, System_Parameters *params) {
	char name[100];
	for (short i = 0; i < 2; i++) {
		fscanf(file, "%s ", name);
		fscanf(file, "%lg ", &params->system[i].bh[0].m);
		fscanf(file, "%lg ", &params->system[i].bh[1].m);
		fscanf(file, "%lg ", &params->system[i].bh[0].chi_Amp);
		fscanf(file, "%lg ", &params->system[i].bh[0].kappa);
		params->system[i].bh[0].kappa *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
		fscanf(file, "%lg ", &params->system[i].bh[0].psi);
		params->system[i].bh[0].psi *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
		fscanf(file, "%lg ", &params->system[i].bh[1].chi_Amp);
		fscanf(file, "%lg ", &params->system[i].bh[1].kappa);
		params->system[i].bh[1].kappa *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
		fscanf(file, "%lg ", &params->system[i].bh[1].psi);
		params->system[i].bh[1].psi *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
		fscanf(file, "%lg ", &params->system[i].dist);
		fscanf(file, "%lg ", &params->system[i].incl);
		params->system[i].incl *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
		fscanf(file, "%lg ", &params->freq_Sampling);
		params->time_Sampling = 1.0 / params->freq_Sampling;
		fscanf(file, "%lg ", &params->freq_Initial);
		fscanf(file, "%lg ", &params->freq_Max);
		fscanf(file, "%s ", params->phase[i]);
		fscanf(file, "%s ", params->spin[i]);
		fscanf(file, "%hd ", &params->amp_Code[i]);
		fscanf(file, "%s\n", params->approx[i]);
		params->system[i].F.declination = params->system[i].F.gmst
				= params->system[i].F.greenwich_Hour_Angle = params->system[i].F.polarization
						= params->system[i].F.right_Ascention = 0.;
		convert_Masses(&params->system[i], FROM_M1M2);
		convert_Spins(&params->system[i], FROM_KAPPA_PSI);
		calc_Antenna_Pattern_For(LH, &params->system[i].F);
	}
}

void multirun(Program_Parameters *program_Parameters, char *file_Name) {
	System_Parameters parameters;
	FILE *file = safely_Open_File_For_Reading(file_Name);
	char temp[FILE_NAME_LENGTH];
	sprintf(temp, "%s", program_Parameters->folder);
	while (!feof(file)) {
		readExactParameters(file, &parameters);
		calc_Matches_For_ParameterPair(program_Parameters, &parameters);
		sprintf(program_Parameters->folder, "%s", temp);
	}
	fclose(file);
}

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
	//binary_System limits_Of_Parameters[2];
	System_Parameters parameters;
	puts("Start!!");
	read_Program_Parameters(&program_Parameters, &parameters, program_Parameters_File_Name);
	multirun(&program_Parameters, parameters_File_Name);
	//read_Parameters(limits_Of_Parameters, parameters_File_Name);
	//proba(&program_Parameters, &parameters, limits_Of_Parameters);
	//run_Algorithm(&program_Parameters, &parameters, limits_Of_Parameters);
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
