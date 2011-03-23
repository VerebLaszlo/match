/**
 * @file main_SQT-ST.c
 *
 * @date Mar 18, 2011
 * @author vereb
 */

#include <stdlib.h>
#include <time.h>
#include "io_handler.h"

#include <lal/LALSQTPNWaveformInterface.h>
#include <lal/LALNoiseModelsInspiral.h>
#include <lal/GenerateInspiral.h>
#include <lal/LALInspiral.h>

int lalDebugLevel = 0;

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

void generate_Parameters1(System_Parameters *parameters, binary_System *limits);

short run_For_Time(Program_Parameters *prog, System_Parameters *parameters);

void initLALParameters1(LALParameters *lalparams, System_Parameters *parameters);

void create_Signal_Struct1(signalStruct *signal, long size);

void set_hphc(short index, long elem, signalStruct *sig, LALParameters *lal,
		System_Parameters *params);

void destroy_Signal_Struct1(signalStruct *signal);

void XLALSQTPNDestroyCoherentGW1(CoherentGW *wave);

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
	generate_Parameters1(&parameters, limits_Of_Parameters);
	run_For_Time(&program_Parameters, &parameters);
	puts("Done!!!");
	return EXIT_SUCCESS;
}

void generate_Parameters1(System_Parameters *parameters, binary_System *limits) {
	assert(parameters);
	assert(limits);
	gen_Parameters(&parameters->system[0], &limits[0], &limits[1], ETAM, KAPPA_PSI);
	memcpy(&parameters->system[1], &parameters->system[0], sizeof(binary_System));
}

short run_For_Time(Program_Parameters *prog, System_Parameters *parameters) {
	assert(prog);
	assert(parameters);
	static LALParameters lalparams;
	signalStruct sig;
	initLALParameters1(&lalparams, parameters);
	for (short i = 0; i < 2; i++) {
		memset(&lalparams.status, 0, sizeof(LALStatus));
		LALGenerateInspiral(&lalparams.status, &lalparams.waveform[i], &lalparams.injParams[i],
				&lalparams.ppnParams);
		if (lalparams.status.statusCode) {
			fprintf(stderr, "%d: LALSQTPNWaveformTest: error generating waveform\n", i);
			XLALSQTPNDestroyCoherentGW1(&lalparams.waveform[0]);
			XLALSQTPNDestroyCoherentGW1(&lalparams.waveform[1]);
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
			set_hphc(i, j, &sig, &lalparams, parameters);
		}
	}
	for (long i = 0; i < parameters->max_Length; i++) {
		double x = fabs(sig.signal[0][i] - sig.signal[2][i] + sig.signal[1][i] - sig.signal[3][i]);
		if (x > 5.0 * 1e-26) {
			//printf("%ld: % 20.15lg\n", i, x);
		}
	}
	write_Wave_To_File1(prog, parameters, &sig, 0);
	XLALSQTPNDestroyCoherentGW1(&lalparams.waveform[0]);
	XLALSQTPNDestroyCoherentGW1(&lalparams.waveform[1]);
	destroy_Signal_Struct1(&sig);
	return 0;
}

void initLALParameters1(LALParameters *lalparams, System_Parameters *parameters) {
	assert(lalparams);
	assert(parameters);
	memset(&lalparams->waveform, 0, 2 * sizeof(CoherentGW));
	memset(&lalparams->injParams, 0, 2 * sizeof(SimInspiralTable));
	memset(&lalparams->ppnParams, 0, sizeof(PPNParamStruc));
	memset(&lalparams->waveform, 0, 2 * sizeof(CoherentGW));
	lalparams->ppnParams.deltaT = 1. / parameters->freq_Sampling;
	parameters->freq_Min = 40.;
	for (short i = 0; i < 2; i++) {
		lalparams->injParams[i].mass1 = parameters->system[i].bh[0].m;
		lalparams->injParams[i].mass2 = parameters->system[i].bh[1].m;
		lalparams->injParams[i].spin1x = parameters->system[i].bh[0].chi[0];
		lalparams->injParams[i].spin1y = parameters->system[i].bh[0].chi[1];
		lalparams->injParams[i].spin1z = parameters->system[i].bh[0].chi[2];
		lalparams->injParams[i].spin2x = parameters->system[i].bh[1].chi[0];
		lalparams->injParams[i].spin2y = parameters->system[i].bh[1].chi[1];
		lalparams->injParams[i].spin2z = parameters->system[i].bh[1].chi[2];
		lalparams->injParams[i].inclination = parameters->system[i].incl;
		lalparams->injParams[i].f_lower = parameters->freq_Min;
		lalparams->injParams[i].distance = parameters->system[i].dist;
		lalparams->injParams[i].coa_phase = parameters->system[i].coaPhase = 0.;
		lalparams->injParams[i].f_lower = parameters->freq_Initial;
		lalparams->ppnParams.deltaT = 1. / parameters->freq_Sampling;
		lalparams->injParams[i].amp_order = parameters->amp_Code[i];
		snprintf(lalparams->injParams[i].waveform, LIGOMETA_WAVEFORM_MAX * sizeof(CHAR), "%s%s%s",
				parameters->approx[i], parameters->phase[i], parameters->spin[i]);
	}
}

void create_Signal_Struct1(signalStruct *signal, long size) {
	assert(size>0);
	signal->size = size;
	short i;
	for (i = 0; i < NUM_OF_SIGNALS; i++) {
		signal->signal[i] = malloc(signal->size * sizeof(double));
		memset(signal->signal[i], 0, signal->size * sizeof(double));
		signal->product_Signal[i] = NULL;
		signal->csignal[i] = NULL;
		signal->plan[i] = NULL;
	}
	signal->psd = NULL;
}

void set_hphc(short index, long elem, signalStruct *sig, LALParameters *lal,
		System_Parameters *params) {
	double a1, a2, phi, shift;
	if (!strcmp(params->approx[index], "SpinQuadTaylor")) {
		sig->signal[2 * index][elem] = lal->waveform[index].h->data->data[2 * elem];
		sig->signal[2 * index + 1][elem] = lal->waveform[index].h->data->data[2 * elem + 1];
	} else if (!strcmp(params->approx[index], "SpinTaylor")) {
		a1 = lal->waveform[index].a->data->data[2 * elem];
		a2 = lal->waveform[index].a->data->data[2 * elem + 1];
		phi = lal->waveform[index].phi->data->data[elem] - lal->waveform[index].phi->data->data[0];
		shift = lal->waveform[index].shift->data->data[elem];
		sig->signal[2 * index][elem] = a1 * cos(shift) * cos(phi) - a2 * sin(shift) * sin(phi);
		sig->signal[2 * index + 1][elem] = a1 * sin(shift) * cos(phi) + a2 * cos(shift) * sin(phi);
	}
}

void destroy_Signal_Struct1(signalStruct *signal) {
	assert(signal);
	short i;
	for (i = 0; i < NUM_OF_SIGNALS; i++) {
		if (signal->signal[i]) {
			free(signal->signal[i]);
		}
		if (signal->product_Signal[i]) {
			free(signal->product_Signal[i]);
		}
		if (signal->csignal[i]) {
			free(signal->csignal[i]);
		}
	}
	if (signal->psd) {
		free(signal->psd);
	}
}

void XLALSQTPNDestroyCoherentGW1(CoherentGW *wave) {
	if (wave->h) {
		if (wave->h->data) {
			XLALDestroyREAL4VectorSequence(wave->h->data);
		}
		XLALFree(wave->h);
		wave->h = NULL;
	}
	if (wave->a) {
		if (wave->a->data) {
			XLALDestroyREAL4VectorSequence(wave->a->data);
		}
		XLALFree(wave->a);
		wave->a = NULL;
	}
	if (wave->f) {
		if (wave->f->data) {
			XLALDestroyREAL4Vector(wave->f->data);
		}
		XLALFree(wave->f);
		wave->f = NULL;
	}
	if (wave->phi) {
		if (wave->phi->data) {
			XLALDestroyREAL8Vector(wave->phi->data);
		}
		XLALFree(wave->phi);
		wave->phi = NULL;
	}
	if (wave->shift) {
		if (wave->shift->data) {
			XLALDestroyREAL4Vector(wave->shift->data);
		}
		XLALFree(wave->shift);
		wave->shift = NULL;
	}
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
