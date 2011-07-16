/*
 * @file match.c
 * @author László Veréb
 * @date 2010.04.08.
 */

#include "match.h"

/**
 * Formulae given.
 * \f{gather*}{
 * 	\inProd{h_1}{h_2}=4\Re\int_{f_{min}}^{f_{max}}\frac{\tilde{h}_1(f)\tilde{h}_2^*(f)}{S_h(f)}df
 * \f}
 * @param left
 * @param right
 * @param norm
 * @param min_Index
 * @param max_Index
 * @return
 */
static double inner_Product(fftw_complex left[], fftw_complex right[], double norm[],
		long min_Index, long max_Index) {
	assert(0<min_Index && min_Index< max_Index);
	double scalar = 0.;
	for (long i = min_Index; i < max_Index; i++) {
		assert(norm[i]);
		scalar += (left[i][0] * right[i][0] + left[i][1] * right[i][1]) / norm[i];
	}
	return 4.0 * scalar;
}

/**
 * Formulae given.
 * \f{gather*}{
 * 	\tilde{e}=\frac{\tilde{h}}{\sqrt{\inProd{\tilde{h}}{\tilde{h}}}}
 * \f}
 * @param in
 * @param out
 * @param min_Index
 * @param max_Index
 */
static void normalise(signalStruct *in, long min_Index, long max_Index, signalStruct *out) {
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

/**
 * Needs formulae.
 * \f{gather*}{
 * 	\tilde{e}_\bot=\tilde{e}_\times-\tilde{e}_+\frac{\inProd{\tilde{e}_+}{\tilde{e}_\times}}{\inProd{\tilde{e}_+}{\tilde{e}_+}}
 * \f}
 * @param in
 * @param out
 * @param min_Index
 * @param max_Index
 */
static void orthogonise(signalStruct *in, long min_Index, long max_Index, signalStruct *out) {
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

/**
 * Needs formulae.
 * @param in
 * @param out
 * @param min_Index
 * @param max_Index
 */
static void orthonormalise(signalStruct *in, long min_Index, long max_Index, signalStruct *out) {
	assert(in);
	assert(out);
	assert(in->size >0 && out->size == in->size);
	assert(0<min_Index && min_Index < max_Index);
	normalise(in, min_Index, max_Index, in);
	orthogonise(in, min_Index, max_Index, out);
}

/**
 * Needs formulae.
 * @param left
 * @param right
 * @param norm
 * @param min_Index
 * @param max_Index
 * @param out
 */
static void cross_Product(fftw_complex left[], fftw_complex right[], double norm[], long min_Index,
		long max_Index, fftw_complex out[]) {
	assert(0<min_Index && min_Index<max_Index);
	for (long i = min_Index; i < max_Index; i++) {
		assert(norm[i]);
		out[i][0] = 4.0 * (left[i][0] * right[i][0] + left[i][1] * right[i][1]) / norm[i];
		out[i][1] = 4.0 * (left[i][1] * right[i][0] - left[i][0] * right[i][1]) / norm[i];
	}
}

/**
 * Needs formulae.
 * @param in
 * @param typ
 * @param best
 * @param minimax
 */
static void calc_Timemaximised_Matches(signalStruct *in, double *typ, double *best, double *minimax) {
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
	calc_Timemaximised_Matches(in, typ, best, minimax);
}

char apr[2][FILENAME_MAX];

void setSignal_From_A1A2(short i, signalStruct *sig, LALParameters *lal, double F[]) {
	REAL8 a1, a2, phi, shift;
	for (long j = 0; j < sig->length[i]; j++) {
		a1 = lal->waveform[i].a->data->data[2 * j];
		a2 = lal->waveform[i].a->data->data[2 * j + 1];
		phi = lal->waveform[i].phi->data->data[j];
		if (strstr(apr[i], "SpinQuadTaylor")) {
			//shift = lal->waveform[i].shift->data->data[j];
			shift = 0.0;
			//phi-=M_PI_2;
		} else if (strstr(apr[i], "SpinTaylorFrameless")) {
			shift = 0.0;
		} else {
			shift = lal->waveform[i].shift->data->data[j];
		}
		//printf("%ld, %ld: %lg %lg %lg %lg\n",sig->length[i],j, a1, a2, phi, shift);
		/*
		 sig->signal[2 * i][j] = lal->waveform[i].a->data->data[2 * j] * M_SQRT2 / 2.0;
		 sig->signal[2 * i + 1][j] = lal->waveform[i].a->data->data[2 * j + 1] * M_SQRT2 / 2.0;
		 sig->signal[RESPONSE1 + i][j] = sig->signal[2 * i][j] * F[0] + sig->signal[2 * i + 1][j]
		 * F[1];
		 */
		sig->signal[2 * i][j] = a1 * cos(shift) * cos(phi) - a2 * sin(shift) * sin(phi);
		sig->signal[2 * i + 1][j] = a1 * sin(shift) * cos(phi) + a2 * cos(shift) * sin(phi);
		sig->signal[RESPONSE1 + i][j] = sig->signal[2 * i][j] * F[0] + sig->signal[2 * i + 1][j]
				* F[1];
	}
}

void setSignal_From_HPHC(short i, signalStruct *sig, LALParameters *lal, double F[]) {
	for (unsigned long j = 0; j < lal->waveform[i].f->data->length; j++) {
		sig->signal[2 * i][j] = lal->waveform[i].h->data->data[2 * j];
		sig->signal[2 * i + 1][j] = lal->waveform[i].h->data->data[2 * j + 1];
		sig->signal[RESPONSE1 + i][j] = sig->signal[2 * i][j] * F[0] + sig->signal[2 * i + 1][j]
				* F[1];
	}
}
