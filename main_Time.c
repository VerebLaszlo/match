/**
 * @file main_Time.c
 *
 * @date Mar 17, 2011
 * @author vereb
 */

#include <stdlib.h>
#include <time.h>

#include "io_handler.h"

int lalDebugLevel = 0;

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

time_t start, end;

void set_hphc(short index, long elem, signalStruct *sig, LALParameters *lal,
		System_Parameters *params) {
	double a1, a2, phi, shift;
	for (short i = 0; i < 2; i++) {
		if (!strcmp(params->approx[i], "SpinQuadTaylor")) {
			sig->signal[2 * index][elem] = lal->waveform[index].h->data->data[2 * elem];
			sig->signal[2 * index + 1][elem] = lal->waveform[index].h->data->data[2 * elem + 1];
		} else if (!strcmp(params->approx[i], "SpinTaylor")) {
			a1 = lal->waveform[index].a->data->data[2 * elem];
			a2 = lal->waveform[index].a->data->data[2 * elem + 1];
			phi = lal->waveform[index].phi->data->data[elem]-lal->waveform[index].phi->data->data[0];
			shift = lal->waveform[index].shift->data->data[elem];
			sig->signal[2 * index][elem] = a1 * cos(shift) * cos(phi) - a2 * sin(shift) * sin(phi);
			sig->signal[2 * index + 1][elem] = a1 * sin(shift) * cos(phi) + a2 * cos(shift) * sin(
					phi);
		}
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
		//		if (signal->plan[i]) {
		//			fftw_destroy_plan(signal->plan[i]);
		//		}
	}
	if (signal->psd) {
		free(signal->psd);
	}
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

void generate_Parameters1(System_Parameters *parameters, binary_System *limits) {
	assert(parameters);
	assert(limits);
	gen_Parameters(&parameters->system[0], &limits[0], &limits[1], ETAM, KAPPA_PSI);
	memcpy(&parameters->system[1], &parameters->system[0], sizeof(binary_System));
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
		}/*
		 if (parameters->shorter) {
		 for (; i < parameters->max_Length; i++) {
		 fprintf(file, "%*.*lg %% ", prog->width_Of_Number_To_Plot, prog->precision_To_Plot,
		 (double)i * parameters->time_Sampling);
		 fprintf(file, text, sig->signal[H1P][i], sig->signal[H1C][i], sig->signal[H1P][i]
		 * parameters->system[0].F.antenna_Beam_Pattern[0] + sig->signal[H1C][i]
		 * parameters->system[0].F.antenna_Beam_Pattern[1]);
		 fprintf(file, "%*s %% %*s %% %*s %% ", prog->width_Of_Number_To_Plot, "",
		 prog->width_Of_Number_To_Plot, "", prog->width_Of_Number_To_Plot, "");
		 fprintf(file, "\n");
		 }
		 } else {
		 fprintf(file, "%*.*lg %% ", prog->width_Of_Number_To_Plot, prog->precision_To_Plot,
		 (double)i * parameters->time_Sampling);
		 fprintf(file, "%*s %*s %*s %%%% ", prog->width_Of_Number_To_Plot, "",
		 prog->width_Of_Number_To_Plot, "", prog->width_Of_Number_To_Plot, "");
		 fprintf(file, text, sig->signal[H2P][i], sig->signal[H2C][i], sig->signal[H2P][i]
		 * parameters->system[0].F.antenna_Beam_Pattern[0] + sig->signal[H2C][i]
		 * parameters->system[0].F.antenna_Beam_Pattern[1]);
		 fprintf(file, "\n");
		 }*/
		fclose(file);
	} else {
		printf("Can not open file: %s, terminating!!!\n", file_Name);
		fflush(stdout);
		abort();
	}
}

short run_For_Time(Program_Parameters *prog, System_Parameters *parameters, long index) {
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
			set_hphc(i, j, &sig, &lalparams, parameters);
			//sig.signal[2 * i][j] = lalparams.waveform[i].h->data->data[2 * j];
			//sig.signal[2 * i + 1][j] = lalparams.waveform[i].h->data->data[2 * j + 1];
		}
	}
	write_Wave_To_File1(prog, parameters, &sig, index);
	XLALSQTPNDestroyCoherentGW1(&lalparams.waveform[0]);
	XLALSQTPNDestroyCoherentGW1(&lalparams.waveform[1]);
	destroy_Signal_Struct1(&sig);
	return 0;
}

void proba1(Program_Parameters *program_Parameters, System_Parameters *parameters,
		binary_System *limits) {
	assert(program_Parameters);
	assert(parameters);
	assert(limits);
	assert(program_Parameters->number_Of_Runs >= 0);
	char temp[FILE_NAME_LENGTH];
	srand(86);
	sprintf(temp, "%s", program_Parameters->folder);
	time(&start);
	for (long i = 0; i < program_Parameters->number_Of_Runs; i++) {
		generate_Parameters1(parameters, limits);
		run_For_Time(program_Parameters, parameters, i);
		if ((i + 1) % 10 == 0) {
			time(&end);
			printf("%ld %lg\n", i + 1, difftime(end, start));
		}
	}
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
	//LALCheckMemoryLeaks();
	puts("Done!!!");
	return EXIT_SUCCESS;
}
