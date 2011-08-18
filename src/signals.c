/*
 * datatypes.c
 *
 *  Created on: Jul 16, 2011
 *      Author: vereb
 */

#include "signals.h"

void create_Signal_Struct(SignalStruct *signal, long size) {
	assert(size>0);
	signal->size = size;
	short i;
	for (i = 0; i < NUMBER_OF_SIGNALS_COMPONENTS; i++) {
		signal->componentsInTime[i] = fftw_malloc(signal->size * sizeof(double));
		memset(signal->componentsInTime[i], 0, signal->size * sizeof(double));
		signal->product[i] = fftw_malloc(signal->size * sizeof(double));
		memset(signal->product[i], 0, signal->size * sizeof(double));
		signal->componentsInFrequency[i] = fftw_malloc(signal->size * sizeof(fftw_complex));
		memset(signal->componentsInFrequency[i], 0, signal->size * sizeof(fftw_complex));
		signal->plan[i] = fftw_plan_dft_r2c_1d(signal->size, signal->componentsInTime[i],
			signal->componentsInFrequency[i], FFTW_ESTIMATE);
	}
	for (; i < NUMBER_OF_SIGNALS; i++) {
		signal->inTime[i] = fftw_malloc(signal->size * sizeof(double));
		memset(signal->inTime[i], 0, signal->size * sizeof(double));
	}
	signal->powerSpectrumDensity = fftw_malloc(signal->size * sizeof(double));
	memset(signal->powerSpectrumDensity, 0, signal->size * sizeof(double));
}

void create_Signal_Struct1(SignalStruct *signal, long size) {
	assert(size>0);
	signal->size = size;
	short i;
	for (i = 0; i < NUMBER_OF_SIGNALS_COMPONENTS; i++) {
		signal->componentsInTime[i] = malloc(signal->size * sizeof(double));
		memset(signal->componentsInTime[i], 0, signal->size * sizeof(double));
		signal->product[i] = NULL;
		signal->componentsInFrequency[i] = NULL;
		signal->plan[i] = NULL;
	}
	for (; i < NUMBER_OF_SIGNALS; i++) {
		signal->inTime[i] = malloc(signal->size * sizeof(double));
		memset(signal->inTime[i], 0, signal->size * sizeof(double));
	}
	signal->powerSpectrumDensity = NULL;
}

void destroy_Signal_Struct(SignalStruct *signal) {
	assert(signal);
	short i;
	for (i = 0; i < NUMBER_OF_SIGNALS_COMPONENTS; i++) {
		if (signal->componentsInTime[i]) {
			fftw_free(signal->componentsInTime[i]);
		}
		if (signal->product[i]) {
			fftw_free(signal->product[i]);
		}
		if (signal->componentsInFrequency[i]) {
			fftw_free(signal->componentsInFrequency[i]);
		}
		if (signal->plan[i]) {
			fftw_destroy_plan(signal->plan[i]);
		}
	}
	for (; i < NUMBER_OF_SIGNALS; i++) {
		fftw_free(signal->inTime[i]);
	}
	if (signal->powerSpectrumDensity) {
		fftw_free(signal->powerSpectrumDensity);
	}
}

void destroy_Signal_Struct1(SignalStruct *signal) {
	assert(signal);
	short i;
	for (i = 0; i < NUMBER_OF_SIGNALS_COMPONENTS; i++) {
		if (signal->componentsInTime[i]) {
			free(signal->componentsInTime[i]);
		}
		if (signal->product[i]) {
			free(signal->product[i]);
		}
		if (signal->componentsInFrequency[i]) {
			free(signal->componentsInFrequency[i]);
		}
	}
	for (; i < NUMBER_OF_SIGNALS; i++) {
		free(signal->inTime[i]);
	}
	if (signal->powerSpectrumDensity) {
		free(signal->powerSpectrumDensity);
	}
}

void print_Two_Signals(FILE*file, SignalStruct *sig, double dt, short width, short precision) {
	short shorter = sig->length[0] < sig->length[1] ? 0 : 1;
	static char temp[FILENAME_MAX];
	static char format[FILENAME_MAX];
	static char text[FILENAME_MAX];
	static char text1[FILENAME_MAX];
	static char textformat[FILENAME_MAX];
	sprintf(temp, "%%%d.%dlg %%%% ", width, precision);
	sprintf(text, "%%%ds %%%% ", width);
	sprintf(format, "%s %s", temp, temp);
	sprintf(text1, "%%-%ds ", width);
	sprintf(textformat, "%s %s %s", text1, text1, text1);
	fprintf(file, textformat, "#time", " h1", "  h2");
	fprintf(file, "\n");
	long i;
	for (i = 0; i < sig->length[shorter]; i++) {
		fprintf(file, temp, (double) i * dt);
		fprintf(file, format, sig->inTime[H1][i], sig->inTime[H2][i]);
		fprintf(file, "\n");
	}
	if (shorter) {
		for (; i < sig->size; i++) {
			fprintf(file, temp, (double) i * dt);
			fprintf(file, temp, sig->inTime[H1][i]);
			fprintf(file, "\n");
		}
	} else {
		for (; i < sig->size; i++) {
			fprintf(file, temp, (double) i * dt);
			fprintf(file, text, "");
			fprintf(file, temp, sig->inTime[H2][i]);
			fprintf(file, "\n");
		}
	}
}

void print_Two_Signals_And_Difference(FILE*file, SignalStruct *sig, double dt, short width,
	short precision) {
	short shorter = sig->length[0] < sig->length[1] ? 0 : 1;
	static char temp[FILENAME_MAX];
	static char format[FILENAME_MAX];
	static char text[FILENAME_MAX];
	sprintf(temp, "%%%d.%dlg %%%% ", width, precision);
	sprintf(text, "%%%ds %%", width);
	sprintf(format, "%s%s%s", temp, temp, temp);
	puts(format);
	long i;
	for (i = 0; i < sig->length[shorter]; i++) {
		fprintf(file, temp, (double) i * dt);
		fprintf(file, format, sig->inTime[H1][i], sig->inTime[H2][i],
			sig->inTime[H1][i] - sig->inTime[H2][i]);
		fprintf(file, "\n");
	}
	if (shorter) {
		for (; i < sig->size; i++) {
			fprintf(file, temp, (double) i * dt);
			fprintf(file, temp, sig->inTime[H1][i]);
			fprintf(file, text, "");
			fprintf(file, temp, sig->inTime[H1][i]);
			fprintf(file, "\n");
		}
	} else {
		for (; i < sig->size; i++) {
			fprintf(file, temp, (double) i * dt);
			fprintf(file, text, width, "");
			fprintf(file, temp, sig->inTime[H2][i]);
			fprintf(file, temp, -sig->inTime[H2][i]);
			fprintf(file, "\n");
		}
	}
}

void print_Two_Signals_With_HPHC(FILE*file, SignalStruct *sig, double dt, short width,
	short precision) {
	short shorter = sig->length[0] < sig->length[1] ? 0 : 1;
	static char temp[FILENAME_MAX];
	static char format[FILENAME_MAX];
	static char text[FILENAME_MAX];
	static char textformat[FILENAME_MAX];
	sprintf(temp, "%%%d.%dlg %%%%", width, precision);
	sprintf(text, "%%%ds %%%%", width);
	sprintf(format, " %s %s %s", temp, temp, temp);
	sprintf(textformat, " %s %s %s", text, text, text);
	long i;
	for (i = 0; i < sig->length[shorter]; i++) {
		fprintf(file, temp, (double) i * dt);
		fprintf(file, format, sig->inTime[H1][i], sig->componentsInTime[H1P][i],
			sig->componentsInTime[H1C][i]);
		fprintf(file, format, sig->inTime[H2][i], sig->componentsInTime[H2P][i],
			sig->componentsInTime[H2C][i]);
		fprintf(file, "\n");
	}
	if (shorter) {
		for (; i < sig->size; i++) {
			fprintf(file, temp, (double) i * dt);
			fprintf(file, format, sig->inTime[H1][i], sig->componentsInTime[H1P][i],
				sig->componentsInTime[H1C][i]);
			fprintf(file, "\n");
		}
	} else {
		for (; i < sig->size; i++) {
			fprintf(file, temp, (double) i * dt);
			fprintf(file, textformat, "", "", "");
			fprintf(file, format, sig->inTime[H2][i], sig->componentsInTime[H2P][i],
				sig->componentsInTime[H2C][i]);
			fprintf(file, "\n");
		}
	}
}

void calculate_H_From_HPHC(SignalStruct *signal, double *antennaFunction) {
	for (short i = 0; i < 2; i++) {
		for (long j = 0; j < signal->length[i]; j++) {
			signal->inTime[H1 + i][j] = signal->inTime[H1P + 2 * i][j]
				* antennaFunction[H1P + 2 * i]
				+ signal->inTime[H1C + 2 * i][j] * antennaFunction[H1C + 2 * i];
		}
	}
}
