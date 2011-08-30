/**
 * @file signals.c
 *
 * @date 2011.08.18.
 * @author László Veréb
 */

#include "signals.h"

void createSignal(SignalStruct *signal, size_t size) {
	assert(signal);
	assert(size);
	memset(signal, 0, sizeof(SignalStruct));
	signal->size = size;
	for (ushort i = 0; i < NUMBER_OF_SIGNALS; i++) {
		signal->inTime[i] = malloc(signal->size * sizeof(double));
		memset(signal->inTime[i], 0, signal->size * sizeof(double));
	}
	for (ushort i = 0; i < NUMBER_OF_SIGNALS_COMPONENTS; i++) {
		signal->componentsInTime[i] = malloc(signal->size * sizeof(double));
		memset(signal->componentsInTime[i], 0, signal->size * sizeof(double));
	}
}

void destroySignal(SignalStruct *signal) {
	assert(signal);
	signal->size = 0;
	for (ushort i = 0; i < NUMBER_OF_SIGNALS; i++) {
		if (signal->inTime[i]) {
			free(signal->inTime[i]);
		}
	}
	for (ushort i = 0; i < NUMBER_OF_SIGNALS_COMPONENTS; i++) {
		free(signal->componentsInTime[i]);
	}
}

void printTwoSignals(FILE *file, SignalStruct *signal, OutputFormat *format) {
	assert(file);
	assert(signal);
	assert(format);
	ushort firstShorter = signal->length[0] < signal->length[1] ? 1 : 0;
	ushort number = 3;
	ushort length = number * format->widthWithSeparator;
	char formatString[length];
	setFormatEnd(formatString, number, format);
	size_t i;
	for (i = 0; i < signal->length[firstShorter]; i++) {
		fprintf(file, formatString, (double) i * signal->samplingTime, signal->inTime[H1][i],
			signal->inTime[H2][i]);
	}
	if (firstShorter) {
		sprintf(formatString, "%s %%%c %s %%%c %s\n", format->oneNumber, format->separator,
			format->empty, format->separator, format->oneNumber);
		for (; i < signal->size; i++) {
			fprintf(file, formatString, (double) i * signal->samplingTime, "",
				signal->inTime[H2][i]);
		}
	} else {
		sprintf(formatString, "%s %%%c %s %%%c %s\n", format->oneNumber, format->separator,
			format->oneNumber, format->separator, format->empty);
		for (; i < signal->size; i++) {
			fprintf(file, formatString, (double) i * signal->samplingTime, signal->inTime[H1][i],
				"");
		}
	}
}

void printTwoSignalsAndDifference(FILE *file, SignalStruct *signal, OutputFormat *format) {
	assert(file);
	assert(signal);
	assert(format);
	ushort firstShorter = signal->length[0] < signal->length[1] ? 1 : 0;
	ushort number = 4;
	ushort length = number * format->widthWithSeparator;
	char formatString[length];
	setFormatEnd(formatString, number, format);
	size_t i;
	for (i = 0; i < signal->length[firstShorter]; i++) {
		fprintf(file, formatString, (double) i * signal->samplingTime, signal->inTime[H1][i],
			signal->inTime[H2][i], signal->inTime[H1][i] - signal->inTime[H2][i]);
	}
	if (firstShorter) {
		sprintf(formatString, "%s %%%c %s %%%c %s %%%c %s\n", format->oneNumber, format->separator,
			format->empty, format->separator, format->oneNumber, format->separator,
			format->oneNumber);
		for (; i < signal->size; i++) {
			fprintf(file, formatString, (double) i * signal->samplingTime, "",
				signal->inTime[H2][i], -signal->inTime[H2][i]);
		}
	} else {
		sprintf(formatString, "%s %%%c %s %%%c %s %%%c %s\n", format->oneNumber, format->separator,
			format->oneNumber, format->separator, format->empty, format->separator,
			format->oneNumber);
		for (; i < signal->size; i++) {
			fprintf(file, formatString, (double) i * signal->samplingTime, signal->inTime[H1][i],
				"", signal->inTime[H1][i]);
		}
	}
}

void printTwoSignalsWithHPHC(FILE* file, SignalStruct *signal, OutputFormat *format) {
	assert(file);
	assert(signal);
	assert(format);
	short firstShorter = signal->length[0] < signal->length[1] ? 1 : 0;
	ushort number = 7;
	ushort length = number * format->widthWithSeparator;
	char formatString[length];
	setFormatEnd(formatString, number, format);
	size_t i;
	for (i = 0; i < signal->length[firstShorter]; i++) {
		fprintf(file, formatString, (double) i * signal->samplingTime, signal->inTime[H1][i],
			signal->inTime[H2][i], signal->componentsInTime[H1P][i],
			signal->componentsInTime[H1C][i], signal->componentsInTime[H2P][i],
			signal->componentsInTime[H2C][i]);
	}
	if (firstShorter) {
		sprintf(formatString, "%s %%%c %s %%%c %s %%%c %s %%%c %s %%%c %s %%%c %s\n",
			format->oneNumber, format->separator, format->empty, format->separator,
			format->oneNumber, format->separator, format->empty, format->separator, format->empty,
			format->separator, format->oneNumber, format->separator, format->oneNumber);
		for (; i < signal->size; i++) {
			fprintf(file, formatString, (double) i * signal->samplingTime, "",
				signal->inTime[H2][i], "", "", signal->componentsInTime[H2P][i],
				signal->componentsInTime[H2C][i]);
		}
	} else {
		sprintf(formatString, "%s %%%c %s %%%c %s %%%c %s %%%c %s %%%c %s %%%c %s\n",
			format->oneNumber, format->separator, format->oneNumber, format->separator,
			format->empty, format->separator, format->oneNumber, format->separator,
			format->oneNumber, format->separator, format->empty, format->separator, format->empty);
		for (; i < signal->size; i++) {
			fprintf(file, formatString, (double) i * signal->samplingTime, signal->inTime[H1][i],
				"", signal->componentsInTime[H1P][i], signal->componentsInTime[H1C][i], "", "");
		}
	}
}

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

void calculate_H_From_HPHC(SignalStruct *signal, double *antennaFunction) {
	for (ushort i = 0; i < 2; i++) {
		for (ulong j = 0; j < signal->length[i]; j++) {
			signal->inTime[H1 + i][j] = signal->inTime[H1P + 2 * i][j]
				* antennaFunction[H1P + 2 * i]
				+ signal->inTime[H1C + 2 * i][j] * antennaFunction[H1C + 2 * i];
		}
	}
}
