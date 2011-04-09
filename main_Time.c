/**
 * @file main_Time.c
 *
 * @date Mar 17, 2011
 * @author vereb
 */

#include "match_qmss.h"

int lalDebugLevel = 0;

time_t start, end;

void proba1(Program_Parameters *program_Parameters, System_Parameters *parameters,
		binary_System *limits);

short run_For_Time(Program_Parameters *prog, System_Parameters *parameters, long index) {
	assert(prog);
	assert(parameters);
	LALParameters lalparams;
	signalStruct sig;
	initLALParameters(&lalparams, parameters);
	for (short i = 0; i < 2; i++) {
		memset(&lalparams.status, 0, sizeof(LALStatus));
		LALGenerateInspiral(&lalparams.status, &lalparams.waveform[i], &lalparams.injParams[i],
				&lalparams.ppnParams);
		if (lalparams.status.statusCode) {
			fprintf(stderr, "%d: LALSQTPNWaveformTest: error generating waveform %d\n", i,
					lalparams.status.statusCode);
			XLALSQTPNDestroyCoherentGW(&lalparams.waveform[0]);
			XLALSQTPNDestroyCoherentGW(&lalparams.waveform[1]);
			return 1;
		}
	}
	parameters->max_Length = parameters->min_Length = lalparams.waveform[0].f->data->length
			> lalparams.waveform[1].f->data->length ? lalparams.waveform[0].f->data->length
			: lalparams.waveform[1].f->data->length;
	create_Signal_Struct1(&sig, parameters->max_Length);
	for (short i = 0; i < 2; i++) {
		setSignal_From_A1A2(i, &sig, &lalparams, parameters->system[i].F.antenna_Beam_Pattern);
	}
	if (!index) {
		FILE *file = safely_Open_File_For_Writing("out/fortime.txt");
		print_Two_Signals(file, &sig, parameters->time_Sampling, prog->width_Of_Number_To_Plot,
				prog->precision_To_Plot);
		fclose(file);
	}
	XLALSQTPNDestroyCoherentGW(&lalparams.waveform[0]);
	XLALSQTPNDestroyCoherentGW(&lalparams.waveform[1]);
	destroy_Signal_Struct1(&sig);
	return 0;
}

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
	FILE*file;
	puts("Start!!");
	file = safely_Open_File_For_Reading(program_Parameters_File_Name);
	read_Program_Parameters(file, &program_Parameters);
	fclose(file);
	file = safely_Open_File_For_Reading(parameters_File_Name);
	read_Binary_Parameter_Limits(file, limits_Of_Parameters);
	fclose(file);
	proba1(&program_Parameters, &parameters, limits_Of_Parameters);
	LALCheckMemoryLeaks();
	puts("Done!!!");
	return EXIT_SUCCESS;
}

void proba1(Program_Parameters *program_Parameters, System_Parameters *parameters,
		binary_System *limits) {
	assert(program_Parameters);
	assert(parameters);
	assert(limits);
	assert(program_Parameters->number_Of_Runs >= 0);
	char temp[FILENAME_MAX];
	srand(86);
	sprintf(temp, "%s/%s%s%d.time", program_Parameters->folder, parameters->approx[0],
			parameters->spin[0], parameters->amp_Code[0]);
	FILE *file = safely_Open_File_For_Writing(temp);
	time(&start);
	for (long i = 0; i < program_Parameters->number_Of_Runs; i++) {
		generate_Same_Parameters(parameters, limits);
		run_For_Time(program_Parameters, parameters, i);
		if ((i + 1) % 100 == 0) {
			time(&end);
			fprintf(file, "%ld %lg\n", i + 1, difftime(end, start));
			printf("%ld %lg\n", i + 1, difftime(end, start));
		}
	}
	fclose(file);
}
