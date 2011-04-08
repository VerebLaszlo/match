/*
 * @file match.c
 * @author László Veréb
 * @date 2010.04.08.
 */

#include "match.h"

void calc_Matches(signalStruct *in, long min_Index, long max_Index, double *typ, double *best,
		double *minimax) {
	assert(in);
	assert(in->size);
	assert(0<min_Index && min_Index< max_Index);
	for (short i = 0; i < NUM_OF_SIGNALS; i++) {
		fftw_execute(in->plan[i]);
	}
	orthonormalise(in, min_Index, max_Index, in);
	fftw_complex *product = fftw_malloc(in->size * sizeof(fftw_complex));
	fftw_plan iplan;
	for (short i = 0; i < NUM_OF_SIGNALS; i++) {
		iplan = fftw_plan_dft_c2r_1d(in->size, product, in->product_Signal[i], FFTW_ESTIMATE);
		memset(product, 0, in->size * sizeof(fftw_complex));
		cross_Product(in->csignal[i / 2], in->csignal[i % 2 + 2], in->psd, min_Index, max_Index,
				product);
		fftw_execute(iplan);
		fftw_destroy_plan(iplan);
	}
	fftw_free(product);
	calc_Timemaximised_Matches(in, min_Index, max_Index, typ, best, minimax);
}

void orthonormalise(signalStruct *in, long min_Index, long max_Index, signalStruct *out) {
	assert(in);
	assert(out);
	assert(in->size >0 && out->size == in->size);
	assert(0<min_Index && min_Index < max_Index);
	normalise(in, min_Index, max_Index, in);
	orthogonise(in, min_Index, max_Index, out);
}

void normalise(signalStruct *in, long min_Index, long max_Index, signalStruct *out) {
	assert(in);
	assert(out);
	assert(in->size >0 && out->size == in->size);
	assert(0<min_Index && min_Index < max_Index);
	double normalising_Constant;
	for (short i = 0; i < NUM_OF_SIGNALS; i++) {
		normalising_Constant = sqrt(inner_Product(in->csignal[i], in->csignal[i], in->psd,
				min_Index, max_Index));
		for (long j = 0; j < in->size; j++) {
			assert(normalising_Constant);
			out->csignal[i][j][0] /= normalising_Constant;
			out->csignal[i][j][1] /= normalising_Constant;
		}
	}
}

double inner_Product(fftw_complex left[], fftw_complex right[], double norm[], long min_Index,
		long max_Index) {
	assert(0<min_Index && min_Index< max_Index);
	double scalar = 0.;
	for (long i = min_Index; i < max_Index; i++) {
		assert(norm[i]);
		scalar += (left[i][0] * right[i][0] + left[i][1] * right[i][1]) / norm[i];
	}
	return 4.0 * scalar;
}

void orthogonise(signalStruct *in, long min_Index, long max_Index, signalStruct *out) {
	assert(out);
	assert(in->size >0 && out->size == in->size);
	assert(0<min_Index && min_Index < max_Index);
	double products[NUM_OF_SIGNALS][NUM_OF_SIGNALS];
	for (short i = 0; i < NUM_OF_SIGNALS; i++) {
		for (short j = 0; j < NUM_OF_SIGNALS; j++) {
			products[i][j] = inner_Product(in->csignal[i], in->csignal[j], in->psd, min_Index,
					max_Index);
		}
	}
	for (short i = 0; i < 2; i++) {
		for (long j = 0; j < in->size; j++) {
			out->csignal[2 * i + 1][j][0] = in->csignal[2 * i + 1][j][0] - in->csignal[2 * i][j][0]
					* products[2 * i][2 * i + 1] / products[2 * i][2 * i];
			out->csignal[2 * i + 1][j][1] = in->csignal[2 * i + 1][j][1] - in->csignal[2 * i][j][1]
					* products[2 * i][2 * i + 1] / products[2 * i][2 * i];
		}
	}
}

void cross_Product(fftw_complex left[], fftw_complex right[], double norm[], long min_Index,
		long max_Index, fftw_complex out[]) {
	assert(0<min_Index && min_Index<max_Index);
	for (long i = min_Index; i < max_Index; i++) {
		assert(norm[i]);
		out[i][0] = 4.0 * (left[i][0] * right[i][0] + left[i][1] * right[i][1]) / norm[i];
		out[i][1] = 4.0 * (left[i][1] * right[i][0] - left[i][0] * right[i][1]) / norm[i];
	}
}

void calc_Timemaximised_Matches(signalStruct *in, long min_Index, long max_Index, double *typ,
		double *best, double *minimax) {
	assert(in);
	assert(in->size);
	double A, B, C;
	double match_typ, max_Typ = 0.0;
	double match_best, max_Best = 0.0;
	double match_minimax, max_Minimax = 0.0;
	for (long i = 0; i < in->size; i++) {
		A = SQR(in->product_Signal[HPP][i]) + SQR(in->product_Signal[HPC][i]);
		B = SQR(in->product_Signal[HCP][i]) + SQR(in->product_Signal[HCC][i]);
		C = in->product_Signal[HPP][i] * in->product_Signal[HCP][i] + in->product_Signal[HPC][i]
				* in->product_Signal[HCC][i];
		match_typ = sqrt(A);
		max_Typ = max_Typ > match_typ ? max_Typ : match_typ;
		match_best = sqrt((A + B) / 2. + sqrt(SQR(A - B) / 4. + SQR(C)));
		max_Best = max_Best > match_best ? max_Best : match_best;
		match_minimax = sqrt((A + B) / 2. - sqrt(SQR(A - B) / 4. + SQR(C)));
		max_Minimax = max_Minimax > match_minimax ? max_Minimax : match_minimax;
	}
	*typ = max_Typ / 2.;
	*best = max_Best / 2.;
	*minimax = max_Minimax / 2.;
}

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
	if (signal->psd) {
		free(signal->psd);
	}
}

void print_Two_Signals(FILE*file, signalStruct *sig, double dt, Program_Parameters *prog) {
	short shorter = sig->length[0] < sig->length[1] ? 0 : 1;
	static char temp[FILENAME_MAX];
	static char format[FILENAME_MAX];
	static char text[FILENAME_MAX];
	sprintf(temp, "%%%d.%dlg %%", prog->width_Of_Number_To_Plot, prog->precision_To_Plot);
	sprintf(text, "%%%ds %%", prog->width_Of_Number_To_Plot);
	sprintf(format, "%s %s", temp, temp);
	long i;
	for (i = 0; i < sig->length[shorter]; i++) {
		fprintf(file, temp, (double)i * dt);
		fprintf(file, format, sig->signal[RESPONSE1][i], sig->signal[RESPONSE2]);
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

void print_Two_Signals_And_Difference(FILE*file, signalStruct *sig, double dt,
		Program_Parameters *prog) {
	short shorter = sig->length[0] < sig->length[1] ? 0 : 1;
	static char temp[FILENAME_MAX];
	static char format[FILENAME_MAX];
	static char text[FILENAME_MAX];
	sprintf(temp, "%%%d.%dlg %%", prog->width_Of_Number_To_Plot, prog->precision_To_Plot);
	sprintf(text, "%%%ds %%", prog->width_Of_Number_To_Plot);
	sprintf(format, "%s %s %s", temp, temp, temp);
	long i;
	for (i = 0; i < sig->length[shorter]; i++) {
		fprintf(file, temp, (double)i * dt);
		fprintf(file, format, sig->signal[RESPONSE1][i], sig->signal[RESPONSE2],
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
			fprintf(file, text, prog->width_Of_Number_To_Plot, "");
			fprintf(file, temp, sig->signal[RESPONSE2][i]);
			fprintf(file, temp, -sig->signal[RESPONSE2][i]);
			fprintf(file, "\n");
		}
	}
}

void print_Two_Signals_With_HPHC(FILE*file, signalStruct *sig, double dt, Program_Parameters *prog) {
	short shorter = sig->length[0] < sig->length[1] ? 0 : 1;
	static char temp[FILENAME_MAX];
	static char format[FILENAME_MAX];
	static char text[FILENAME_MAX];
	static char textformat[FILENAME_MAX];
	sprintf(temp, "%%%d.%dlg %%", prog->width_Of_Number_To_Plot, prog->precision_To_Plot);
	sprintf(text, "%%%ds %%", prog->width_Of_Number_To_Plot);
	sprintf(format, "%s %s %s", temp, temp, temp);
	sprintf(textformat, "%s %s %s", text, text, text);
	long i;
	for (i = 0; i < sig->length[shorter]; i++) {
		fprintf(file, temp, (double)i * dt);
		fprintf(file, format, sig->signal[RESPONSE1][i], sig->signal[H1P], sig->signal[H1C]);
		fprintf(file, format, sig->signal[RESPONSE2][i], sig->signal[H2P], sig->signal[H2C]);
		fprintf(file, "\n");
	}
	if (shorter) {
		for (; i < sig->size; i++) {
			fprintf(file, temp, (double)i * dt);
			fprintf(file, format, sig->signal[RESPONSE1][i], sig->signal[H1P], sig->signal[H1C]);
			fprintf(file, "\n");
		}
	} else {
		for (; i < sig->size; i++) {
			fprintf(file, temp, (double)i * dt);
			fprintf(file, textformat, "", "", "");
			fprintf(file, format, sig->signal[RESPONSE2][i], sig->signal[H2P], sig->signal[H2C]);
			fprintf(file, "\n");
		}
	}
}

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

void read_Program_Parameters(Program_Parameters *parameters, System_Parameters *params,
		char *file_Name) {
	assert(parameters);
	assert(params);
	assert(file_Name);
	FILE *file = fopen(file_Name, "r");
	fscanf(file, "%ld\n", &parameters->number_Of_Runs);
	fscanf(file, "%hd\n", &parameters->precision);
	parameters->width_Of_Number = parameters->precision + SPECIAL_CHARACTER_LENGTH;
	fscanf(file, "%hd\n", &parameters->precision_To_Plot);
	parameters->width_Of_Number_To_Plot = parameters->precision_To_Plot + SPECIAL_CHARACTER_LENGTH;
	fscanf(file, "%s\n", parameters->folder);
	fscanf(file, "%lg\n", &params->min_Match);
	fscanf(file, "%lg\n", &params->max_Spin);
	fscanf(file, "%lg\n", &params->spin_Step);
	fscanf(file, "%lg\n", &params->freq_Sampling);
	params->time_Sampling = 1. / params->freq_Sampling;
	fscanf(file, "%lg\n", &params->freq_Initial);
	fscanf(file, "%lg\n", &params->freq_Max);
	fscanf(file, "%lg\n", &params->delta_Length);
	fscanf(file, "%s\n", params->approx[0]);
	fscanf(file, "%s\n", params->phase[0]);
	fscanf(file, "%hd\n", &params->amp_Code[0]);
	fscanf(file, "%s\n", params->spin[0]);
	fscanf(file, "%s\n", params->approx[1]);
	fscanf(file, "%s\n", params->phase[1]);
	fscanf(file, "%hd\n", &params->amp_Code[1]);
	fscanf(file, "%s\n", params->spin[1]);
	fclose(file);
}
