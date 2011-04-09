/**
 * @file main_SQT-ST.c
 *
 * @date Mar 18, 2011
 * @author vereb
 */

#include "match_qmss.h"

int lalDebugLevel = 0;

short run_For_Time(Program_Parameters *prog, System_Parameters *parameters);

int main(int argc, char *argv[]) {
	char program_Parameters_File_Name[FILENAME_MAX];
	char parameters_File_Name[FILENAME_MAX];
	if (argc != 3) {
		puts("\"file for program parameters\" \"file for the parameter limits\"!!!");
		exit(-1);
	}
	sprintf(program_Parameters_File_Name, argv[1]);
	sprintf(parameters_File_Name, argv[2]);
	Program_Parameters program_Parameters;
	binary_System limits_Of_Parameters[2];
	System_Parameters parameters;
	FILE *file;
	puts("Start!!");
	read_Program_Parameters(&program_Parameters, &parameters, program_Parameters_File_Name);
	file = safely_Open_File_For_Reading(parameters_File_Name);
	read_Binary_Parameter_Limits(file, limits_Of_Parameters);
	fclose(file);
	srand(time(NULL));
	srand(86);
	generate_Same_Parameters(&parameters, limits_Of_Parameters);
	run_For_Time(&program_Parameters, &parameters);
	puts("Done!!!");
	return EXIT_SUCCESS;
}

short run_For_Time(Program_Parameters *prog, System_Parameters *parameters) {
	assert(prog);
	assert(parameters);
	static LALParameters lalparams;
	signalStruct sig;
	initLALParameters(&lalparams, parameters);
	for (short i = 0; i < 2; i++) {
		memset(&lalparams.status, 0, sizeof(LALStatus));
		LALGenerateInspiral(&lalparams.status, &lalparams.waveform[i], &lalparams.injParams[i],
				&lalparams.ppnParams);
		if (lalparams.status.statusCode) {
			fprintf(stderr, "%d: LALSQTPNWaveformTest: error generating waveform\n", i);
			XLALSQTPNDestroyCoherentGW(&lalparams.waveform[0]);
			XLALSQTPNDestroyCoherentGW(&lalparams.waveform[1]);
			exit(EXIT_FAILURE);
		}
		puts("X");
	}
	parameters->max_Length = lalparams.waveform[0].f->data->length
			> lalparams.waveform[1].f->data->length ? lalparams.waveform[0].f->data->length
			: lalparams.waveform[1].f->data->length;
	parameters->min_Length = lalparams.waveform[0].f->data->length
			< lalparams.waveform[1].f->data->length ? lalparams.waveform[0].f->data->length
			: lalparams.waveform[1].f->data->length;
	create_Signal_Struct1(&sig, parameters->max_Length);
	for (short i = 0; i < 2; i++) {
		for (unsigned long j = 0; j < lalparams.waveform[i].f->data->length; j++) {
			setSignal_From_A1A2(i, j, &sig, &lalparams,
					parameters->system[i].F.antenna_Beam_Pattern);
		}
	}
	FILE*file = safely_Open_File_For_Writing("out/temp.txt");
	print_Two_Signals(file, &sig, parameters->time_Sampling, prog->width_Of_Number_To_Plot,
			prog->precision_To_Plot);
	XLALSQTPNDestroyCoherentGW(&lalparams.waveform[0]);
	XLALSQTPNDestroyCoherentGW(&lalparams.waveform[1]);
	destroy_Signal_Struct1(&sig);
	return 0;
}
