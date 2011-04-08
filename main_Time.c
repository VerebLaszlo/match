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

void write_Wave_To_File1(Program_Parameters *prog, System_Parameters *parameters,
		signalStruct *sig, short index);

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
		for (long j = 0; j < lalparams.waveform[i].f->data->length; j++) {
			set_Signal_From_A1A2(i, j, &sig, &lalparams);
		}
	}
	if (!index) {
		write_Wave_To_File1(prog, parameters, &sig, index);
	}
	XLALSQTPNDestroyCoherentGW(&lalparams.waveform[0]);
	XLALSQTPNDestroyCoherentGW(&lalparams.waveform[1]);
	destroy_Signal_Struct1(&sig);
	return 0;
}

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
	puts("Start!!");
	read_Program_Parameters(&program_Parameters, &parameters, program_Parameters_File_Name);
	read_Parameters(limits_Of_Parameters, parameters_File_Name);
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
	char temp[2 * FILE_NAME_LENGTH];
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

void write_Wave_To_File1(Program_Parameters *prog, System_Parameters *parameters,
		signalStruct *sig, short index) {
	assert(prog);
	assert(parameters);
	assert(sig);
	char file_Name[FILE_NAME_LENGTH];
	static char temp[FILE_NAME_LENGTH];
	static char text[FILE_NAME_LENGTH];
	sprintf(file_Name, "out/wave%d.txt", index);
	FILE *file = fopen(file_Name, "w");
	if (file) {
		sprintf(temp, "%%%d.%dlg ", prog->width_Of_Number_To_Plot, prog->precision_To_Plot);
		sprintf(text, "%s%s%s", temp, temp, temp);
		fprintf(file, "#");
		fprintf(file, text, parameters->system[0].M, parameters->system[0].bh[0].m
				/ parameters->system[0].bh[1].m, parameters->system[0].eta);
		fprintf(file, text, parameters->system[0].bh[0].chi_Amp, parameters->system[0].bh[0].kappa,
				parameters->system[0].bh[0].psi);
		fprintf(file, text, parameters->system[0].bh[1].chi_Amp, parameters->system[0].bh[1].kappa,
				parameters->system[0].bh[1].psi);
		fprintf(file, "%s %d %s\n", parameters->phase[0], parameters->amp_Code[0],
				parameters->spin[0]);
		fprintf(file, "#");
		fprintf(file, text, parameters->system[1].M, parameters->system[1].bh[0].m
				/ parameters->system[1].bh[1].m, parameters->system[1].eta);
		fprintf(file, text, parameters->system[1].bh[0].chi_Amp, parameters->system[1].bh[0].kappa,
				parameters->system[1].bh[0].psi);
		fprintf(file, text, parameters->system[1].bh[1].chi_Amp, parameters->system[1].bh[1].kappa,
				parameters->system[1].bh[1].psi);
		fprintf(file, "%s %d %s", parameters->phase[1], parameters->amp_Code[1],
				parameters->spin[1]);
		fprintf(file, "\n");
		sprintf(text, "%s %%%% %s %%%% %s %%%% ", temp, temp, temp);
		long i;
		for (i = 0; i < parameters->min_Length; i++) {
			fprintf(file, "%*.*lg %% ", prog->width_Of_Number_To_Plot, prog->precision_To_Plot,
					(double)i * parameters->time_Sampling);
			fprintf(file, text, sig->signal[H1P][i], sig->signal[H1C][i], sig->signal[H1P][i]
					* parameters->system[0].F.antenna_Beam_Pattern[0] + sig->signal[H1C][i]
					* parameters->system[0].F.antenna_Beam_Pattern[1]);
			fprintf(file, text, sig->signal[H2P][i], sig->signal[H2C][i], sig->signal[H2P][i]
					* parameters->system[0].F.antenna_Beam_Pattern[0] + sig->signal[H2C][i]
					* parameters->system[0].F.antenna_Beam_Pattern[1]);
			fprintf(file, "\n");
		}
		fclose(file);
	} else {
		printf("Can not open file: %s, terminating!!!\n", file_Name);
		fflush(stdout);
		abort();
	}
}
