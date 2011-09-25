/**
 * @file match.c
 * @author László Veréb
 * @date 2010.04.08.
 */

#include <math.h>
#include "match.h"
#include "util_math.h"

/**
 * Formulae given.
 * \f[
 * 	\inProd{h_1}{h_2}=4\Re\int_{f_{min}}^{f_{max}}\frac{\tilde{h}_1(f)\tilde{h}_2^*(f)}{S_h(f)}df
 * \f]
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
 * \f[
 * 	\tilde{e}=\frac{\tilde{h}}{\sqrt{\inProd{\tilde{h}}{\tilde{h}}}}
 * \f]
 * @param in
 * @param out
 * @param min_Index
 * @param max_Index
 */
static void normalise(SignalStruct *in, long min_Index, long max_Index, SignalStruct *out) {
	assert(in);
	assert(out);
	assert(in->size >0 && out->size == in->size);
	assert(0<min_Index && min_Index < max_Index);
	double normalising_Constant;
	for (short i = 0; i < NUMBER_OF_SIGNALS_COMPONENTS; i++) {
		normalising_Constant = sqrt(
			inner_Product(in->componentsInFrequency[i], in->componentsInFrequency[i],
				in->powerSpectrumDensity, min_Index, max_Index));
		for (size_t j = 0; j < in->size; j++) {
			assert(normalising_Constant);
			out->componentsInFrequency[i][j][0] /= normalising_Constant;
			out->componentsInFrequency[i][j][1] /= normalising_Constant;
		}
	}
}

/**
 * Needs formulae.
 * \f[
 * 	\tilde{e}_\bot=\tilde{e}_\times-\tilde{e}_+\frac{\inProd{\tilde{e}_+}{\tilde{e}_\times}}{\inProd{\tilde{e}_+}{\tilde{e}_+}}
 * \f]
 * @param in
 * @param out
 * @param min_Index
 * @param max_Index
 */
static void orthogonise(SignalStruct *in, long min_Index, long max_Index, SignalStruct *out) {
	assert(out);
	assert(in->size >0 && out->size == in->size);
	assert(0<min_Index && min_Index < max_Index);
	double products[NUM_OF_SIGNALS][NUM_OF_SIGNALS];
	for (short i = 0; i < NUM_OF_SIGNALS; i++) {
		for (short j = 0; j < NUM_OF_SIGNALS; j++) {
			products[i][j] = inner_Product(in->componentsInFrequency[i],
				in->componentsInFrequency[j], in->powerSpectrumDensity, min_Index, max_Index);
		}
	}
	for (short i = 0; i < 2; i++) {
		for (size_t j = 0; j < in->size; j++) {
			out->componentsInFrequency[2 * i + 1][j][0] = in->componentsInFrequency[2 * i + 1][j][0]
				- in->componentsInFrequency[2 * i][j][0] * products[2 * i][2 * i + 1]
					/ products[2 * i][2 * i];
			out->componentsInFrequency[2 * i + 1][j][1] = in->componentsInFrequency[2 * i + 1][j][1]
				- in->componentsInFrequency[2 * i][j][1] * products[2 * i][2 * i + 1]
					/ products[2 * i][2 * i];
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
static void orthonormalise(SignalStruct *in, long min_Index, long max_Index, SignalStruct *out) {
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
static void calc_Timemaximised_Matches(SignalStruct *in, double *typ, double *best, double *minimax) {
	assert(in);
	assert(in->size);
	double A, B, C;
	double match_typ, max_Typ = 0.0;
	double match_best, max_Best = 0.0;
	double match_minimax, max_Minimax = 0.0;
	for (size_t i = 0; i < in->size; i++) {
		A = square(in->product[HPP][i]) + square(in->product[HPC][i]);
		B = square(in->product[HCP][i]) + square(in->product[HCC][i]);
		C = in->product[HPP][i] * in->product[HCP][i] + in->product[HPC][i] * in->product[HCC][i];
		match_typ = sqrt(A);
		max_Typ = max_Typ > match_typ ? max_Typ : match_typ;
		match_best = sqrt((A + B) / 2. + sqrt(square(A - B) / 4. + square(C)));
		max_Best = max_Best > match_best ? max_Best : match_best;
		match_minimax = sqrt((A + B) / 2. - sqrt(square(A - B) / 4. + square(C)));
		max_Minimax = max_Minimax > match_minimax ? max_Minimax : match_minimax;
	}
	*typ = max_Typ / 2.;
	*best = max_Best / 2.;
	*minimax = max_Minimax / 2.;
}

void calc_Matches(SignalStruct *in, long min_Index, long max_Index, double *typ, double *best,
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
		iplan = fftw_plan_dft_c2r_1d((int) in->size, product, in->product[i], FFTW_ESTIMATE);
		memset(product, 0, in->size * sizeof(fftw_complex));
		cross_Product(in->componentsInFrequency[i / 2], in->componentsInFrequency[i % 2 + 2],
			in->powerSpectrumDensity, min_Index, max_Index, product);
		fftw_execute(iplan);
		fftw_destroy_plan(iplan);
	}
	fftw_free(product);
	calc_Timemaximised_Matches(in, typ, best, minimax);
}
