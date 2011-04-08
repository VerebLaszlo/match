/**
 * @file main_SQT-ST.c
 *
 * @date Mar 18, 2011
 * @author vereb
 */

#include "match_qmss.h"

int lalDebugLevel = 0;

short run_For_Time(Program_Parameters *prog, System_Parameters *parameters);

void write_Wave_To_File1(Program_Parameters *prog, System_Parameters *parameters,
		signalStruct *sig, short index);

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
		for (long j = 0; j < lalparams.waveform[i].f->data->length; j++) {
			set_Signal_From_A1A2(i, j, &sig, &lalparams);
		}
	}
	write_Wave_To_File1(prog, parameters, &sig, 0);
	XLALSQTPNDestroyCoherentGW(&lalparams.waveform[0]);
	XLALSQTPNDestroyCoherentGW(&lalparams.waveform[1]);
	destroy_Signal_Struct1(&sig);
	return 0;
}

void write_Wave_To_File1(Program_Parameters *prog, System_Parameters *parameters,
		signalStruct *sig, short index) {
	assert(prog);
	assert(parameters);
	assert(sig);
	char file_Name[FILE_NAME_LENGTH];
	static char temp[FILE_NAME_LENGTH];
	static char text[FILE_NAME_LENGTH];
	FILE *file;
	for (short i = 0; i < 2; i++) {
		sprintf(file_Name, "out/wave%d.txt", i);
		file = fopen(file_Name, "w");
		sprintf(temp, "%%- %d.%dlg ", prog->width_Of_Number_To_Plot, prog->precision_To_Plot);
		sprintf(text, "%s%s%s", temp, temp, temp);
		fprintf(file, "#");
		fprintf(file, text, parameters->system[i].M, parameters->system[i].bh[0].m
				/ parameters->system[i].bh[1].m, parameters->system[i].eta);
		fprintf(file, text, parameters->system[i].bh[0].chi_Amp, parameters->system[i].bh[0].kappa,
				parameters->system[i].bh[0].psi);
		fprintf(file, text, parameters->system[i].bh[1].chi_Amp, parameters->system[i].bh[1].kappa,
				parameters->system[i].bh[1].psi);
		fprintf(file, "%s %d %s\n", parameters->phase[i], parameters->amp_Code[i],
				parameters->spin[i]);
		fprintf(file, "# ");
		fprintf(file, text, parameters->system[i].F.antenna_Beam_Pattern[0],
				parameters->system[i].F.antenna_Beam_Pattern[1], NAN);
		fprintf(file, "\n");
		sprintf(temp, "%%- %d.%dlg ", prog->width_Of_Number, prog->precision);
		sprintf(text, "%s %%%% %s %%%% %s %%%% ", temp, temp, temp);
		for (long j = 0; j < parameters->min_Length; j++) {
			fprintf(file, "%*.*lg %% ", prog->width_Of_Number_To_Plot, prog->precision_To_Plot,
					(double)j * parameters->time_Sampling);
			fprintf(file, text, sig->signal[2 * i][j], sig->signal[2 * i + 1][j],
					sig->signal[2 * i][j] * parameters->system[i].F.antenna_Beam_Pattern[0]
							+ sig->signal[2 * i + 1][j]
									* parameters->system[i].F.antenna_Beam_Pattern[1]);
			fprintf(file, "\n");
		}
		fclose(file);
	}
}
