/*
 * datatypes.c
 *
 *  Created on: Jul 16, 2011
 *      Author: vereb
 */

#include "datatypes.h"

void create_Signal_Struct(signalStruct *signal, long size) {
	assert(size>0);
	signal->size = size;
	short i;
	for (i = 0; i < NUM_OF_SIGNALS; i++) {
		signal->signal[i] = fftw_malloc(signal->size * sizeof(double));
		memset(signal->signal[i], 0, signal->size * sizeof(double));
		signal->product_Signal[i] = fftw_malloc(signal->size * sizeof(double));
		memset(signal->product_Signal[i], 0, signal->size * sizeof(double));
		signal->csignal[i] = fftw_malloc(signal->size * sizeof(fftw_complex));
		memset(signal->csignal[i], 0, signal->size * sizeof(fftw_complex));
		signal->plan[i] = fftw_plan_dft_r2c_1d(signal->size, signal->signal[i], signal->csignal[i],
				FFTW_ESTIMATE);
	}
	for (; i < NOS_WITH_DETECTOR_RESPONSE; i++) {
		signal->signal[i] = fftw_malloc(signal->size * sizeof(double));
		memset(signal->signal[i], 0, signal->size * sizeof(double));
	}
	signal->psd = fftw_malloc(signal->size * sizeof(double));
	memset(signal->psd, 0, signal->size * sizeof(double));
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
	for (; i < NOS_WITH_DETECTOR_RESPONSE; i++) {
		signal->signal[i] = malloc(signal->size * sizeof(double));
		memset(signal->signal[i], 0, signal->size * sizeof(double));
	}
	signal->psd = NULL;
}

void destroy_Signal_Struct(signalStruct *signal) {
	assert(signal);
	short i;
	for (i = 0; i < NUM_OF_SIGNALS; i++) {
		if (signal->signal[i]) {
			fftw_free(signal->signal[i]);
		}
		if (signal->product_Signal[i]) {
			fftw_free(signal->product_Signal[i]);
		}
		if (signal->csignal[i]) {
			fftw_free(signal->csignal[i]);
		}
		if (signal->plan[i]) {
			fftw_destroy_plan(signal->plan[i]);
		}
	}
	for (; i < NOS_WITH_DETECTOR_RESPONSE; i++) {
		fftw_free(signal->signal[i]);
	}
	if (signal->psd) {
		fftw_free(signal->psd);
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
	for (; i < NOS_WITH_DETECTOR_RESPONSE; i++) {
		free(signal->signal[i]);
	}
	if (signal->psd) {
		free(signal->psd);
	}
}

void print_Two_Signals(FILE*file, signalStruct *sig, double dt, short width, short precision) {
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
		fprintf(file, temp, (double)i * dt);
		fprintf(file, format, sig->signal[RESPONSE1][i], sig->signal[RESPONSE2][i]);
		fprintf(file, "\n");
	}
	if (shorter) {
		for (; i < sig->size; i++) {
			fprintf(file, temp, (double)i * dt);
			fprintf(file, temp, sig->signal[RESPONSE1][i]);
			fprintf(file, "\n");
		}
	} else {
		for (; i < sig->size; i++) {
			fprintf(file, temp, (double)i * dt);
			fprintf(file, text, "");
			fprintf(file, temp, sig->signal[RESPONSE2][i]);
			fprintf(file, "\n");
		}
	}
}

void print_Two_Signals_And_Difference(FILE*file, signalStruct *sig, double dt, short width,
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
		fprintf(file, temp, (double)i * dt);
		fprintf(file, format, sig->signal[RESPONSE1][i], sig->signal[RESPONSE2][i],
				sig->signal[RESPONSE1][i] - sig->signal[RESPONSE2][i]);
		fprintf(file, "\n");
	}
	if (shorter) {
		for (; i < sig->size; i++) {
			fprintf(file, temp, (double)i * dt);
			fprintf(file, temp, sig->signal[RESPONSE1][i]);
			fprintf(file, text, "");
			fprintf(file, temp, sig->signal[RESPONSE1][i]);
			fprintf(file, "\n");
		}
	} else {
		for (; i < sig->size; i++) {
			fprintf(file, temp, (double)i * dt);
			fprintf(file, text, width, "");
			fprintf(file, temp, sig->signal[RESPONSE2][i]);
			fprintf(file, temp, -sig->signal[RESPONSE2][i]);
			fprintf(file, "\n");
		}
	}
}

void print_Two_Signals_With_HPHC(FILE*file, signalStruct *sig, double dt, short width,
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
		fprintf(file, temp, (double)i * dt);
		fprintf(file, format, sig->signal[RESPONSE1][i], sig->signal[H1P][i], sig->signal[H1C][i]);
		fprintf(file, format, sig->signal[RESPONSE2][i], sig->signal[H2P][i], sig->signal[H2C][i]);
		fprintf(file, "\n");
	}
	if (shorter) {
		for (; i < sig->size; i++) {
			fprintf(file, temp, (double)i * dt);
			fprintf(file, format, sig->signal[RESPONSE1][i], sig->signal[H1P][i],
					sig->signal[H1C][i]);
			fprintf(file, "\n");
		}
	} else {
		for (; i < sig->size; i++) {
			fprintf(file, temp, (double)i * dt);
			fprintf(file, textformat, "", "", "");
			fprintf(file, format, sig->signal[RESPONSE2][i], sig->signal[H2P][i],
					sig->signal[H2C][i]);
			fprintf(file, "\n");
		}
	}
}
