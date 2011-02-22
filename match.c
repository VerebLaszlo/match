/**
 * @file match.c
 * @author László Veréb
 * @date 2010.04.09.
 */

#include "match.h"

double pi;///<a

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
		cross_Product(in->csignal[i / 2], in->csignal[i % 2 + 2], in->psd, min_Index, max_Index,
				product);
		iplan = fftw_plan_dft_c2r_1d(in->size, product, in->normalised_Signal[i], FFTW_ESTIMATE);
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
		normalising_Constant = inner_Product(in->csignal[i], in->csignal[i], in->psd, min_Index,
				max_Index);
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

void calc_Timemaximised_Matches(signalStruct *in, long min_Index, long max_Index, double *typ, double *best, double *minimax) {
	assert(in);
	assert(in->size);
	double A, B, C;
	double match_typ, max_Typ;
	double match_best, max_Best;
	double match_minimax, max_Minimax;
	//for (long i = 0; i < in->size; i++) {
		for (long i = min_Index; i < max_Index; i++) {
		A = SQR(in->normalised_Signal[H1P][i]) + SQR(in->normalised_Signal[H1C][i]);
		B = SQR(in->normalised_Signal[H2P][i]) + SQR(in->normalised_Signal[H2C][i]);
		C = in->normalised_Signal[H1P][i] * in->normalised_Signal[H2P][i]
				+ in->normalised_Signal[H1C][i] * in->normalised_Signal[H2C][i];
		match_typ = sqrt(B);
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

// old

int create_Signal_Struct(signalStruct *signal, long size) {
	assert(size>0);
	signal->size = size;
	short i;
	for (i = 0; i < NUM_OF_SIGNALS; i++) {
		signal->signal[i] = fftw_malloc(signal->size * sizeof(double));
		memset(signal->signal[i], 0, signal->size * sizeof(double));
		signal->normalised_Signal[i] = fftw_malloc(signal->size * sizeof(double));
		memset(signal->normalised_Signal[i], 0, signal->size * sizeof(double));
		signal->csignal[i] = fftw_malloc(signal->size * sizeof(fftw_complex));
		memset(signal->csignal[i], 0, signal->size * sizeof(fftw_complex));
		signal->plan[i] = fftw_plan_dft_r2c_1d(signal->size, signal->signal[i], signal->csignal[i],
				FFTW_ESTIMATE);
	}
	signal->psd = fftw_malloc(signal->size * sizeof(double));
	memset(signal->psd, 0, signal->size * sizeof(double));
	return MATCH_SUCCES;
}

void destroy_Signal_Struct(signalStruct *signal) {
	assert(signal);
	short i;
	for (i = 0; i < NUM_OF_SIGNALS; i++) {
		if (signal->signal[i]) {
			fftw_free(signal->signal[i]);
		}
		if (signal->normalised_Signal[i]) {
			fftw_free(signal->normalised_Signal[i]);
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

double inner_Product_Old(fftw_complex left[], fftw_complex right[], double norm[], long minfr,
		long maxfr) {
	assert(minfr>0);
	assert(minfr< maxfr);
	double scalar = 0.;
	long i;
	for (i = minfr; i < maxfr; i++) {
		scalar += ((left[i][0] * right[i][0] + left[i][1] * right[i][1]) / norm[i]);
	}
	return 4. * scalar;
}

void product(fftw_complex out[], fftw_complex left[], fftw_complex right[], double norm[],
		long minfr, long maxfr) {
	assert(minfr>0);
	assert(minfr< maxfr);
	long i;
	for (i = minfr; i < maxfr; i++) {
		if (norm[i] != 0.) {
			out[i][0] = 4. * (left[i][0] * right[i][0] + left[i][1] * right[i][1]) / norm[i];
			out[i][1] = 4. * (left[i][1] * right[i][0] - left[i][0] * right[i][1]) / norm[i];
		} else {
			out[i][0] = out[i][1] = 0.;
		}
	}
}

void normalise_Old(signalStruct *out, signalStruct *in, long minfr, long maxfr) {
	assert(out);
	assert(in);
	assert(minfr>0);
	assert(minfr< maxfr);
	double prod;
	short i;
	long j;
	if (out) {
		// the normalised signals are stored in the out structure
		out->size = in->size;
		for (i = 0; i < NUM_OF_SIGNALS; i++) {
			prod = inner_Product_Old(in->csignal[i], in->csignal[i], in->psd, minfr, maxfr);
			for (j = 0; j < in->size; j++) {
				out->csignal[i][j][0] = in->csignal[i][j][0] / prod;
				out->csignal[i][j][1] = in->csignal[i][j][1] / prod;
			}
		}
	} else {
		// the normalised signals overwrites the original ones
		for (i = 0; i < NUM_OF_SIGNALS; i++) {
			prod = sqrt(inner_Product_Old(in->csignal[i], in->csignal[i], in->psd, minfr, maxfr));
			for (j = 0; j < in->size; j++) {
				in->csignal[i][j][0] /= prod;
				in->csignal[i][j][1] /= prod;
			}
		}
	}
}

void orthogonise_Old(signalStruct *out, signalStruct *in, long minfr, long maxfr) {
	assert(out);
	assert(in);
	assert(minfr>0);
	assert(minfr< maxfr);
	double prod[NUM_OF_SIGNALS][NUM_OF_SIGNALS];
	short i, j;
	long k;
	if (out) {
		// the orthonormalised signals are stored in the out structure
		out->size = in->size;
		double prod[NUM_OF_SIGNALS][NUM_OF_SIGNALS];
		for (i = 0; i < NUM_OF_SIGNALS; i++) {
			for (j = 0; j < NUM_OF_SIGNALS; j++) {
				prod[i][j] = inner_Product_Old(in->csignal[i], in->csignal[j], in->psd, minfr,
						maxfr);
			}
		}
		for (i = 0; i < 2; i++) {
			for (k = 0; k < in->size; k++) {
				out->csignal[2 * i + 1][k][0] = in->csignal[2 * i + 1][k][0]
						- in->csignal[2 * i][k][0] * prod[2 * i][2 * i + 1] / prod[2 * i][2 * i];
				out->csignal[2 * i + 1][k][1] = in->csignal[2 * i + 1][k][1]
						- in->csignal[2 * i][k][1] * prod[2 * i][2 * i + 1] / prod[2 * i][2 * i];
			}
		}
	} else {
		// the orthonormalised signals overwrites the original ones
		for (i = 0; i < NUM_OF_SIGNALS; i++) {
			for (j = 0; j < NUM_OF_SIGNALS; j++) {
				prod[i][j] = inner_Product_Old(in->csignal[i], in->csignal[j], in->psd, minfr,
						maxfr);
			}
		}
		for (i = 0; i < 2; i++) {
			for (k = 0; k < in->size; k++) {
				in->csignal[2 * i + 1][k][0] -= in->csignal[2 * i][k][0] * prod[2 * i][2 * i + 1]
						/ prod[2 * i][2 * i];
				in->csignal[2 * i + 1][k][1] -= in->csignal[2 * i][k][1] * prod[2 * i][2 * i + 1]
						/ prod[2 * i][2 * i];
			}
		}
	}
}

void orthonormalise_Old(signalStruct *out, signalStruct *in, long minfr, long maxfr) {
	assert(in);
	assert(minfr>0);
	assert(minfr< maxfr);
	double prod[NUM_OF_SIGNALS][NUM_OF_SIGNALS];
	short i, j;
	long k;
	if (out) {
		// the orthonormalised signals are stored in the out structure
		out->size = in->size;
		for (i = 0; i < NUM_OF_SIGNALS; i++) {
			out->size = in->size;
			prod[i][i] = sqrt(inner_Product_Old(in->csignal[i], in->csignal[i], in->psd, minfr,
					maxfr));
			for (k = 0; k < in->size; k++) {
				out->csignal[i][k][0] = in->csignal[i][k][0] / prod[i][i];
				out->csignal[i][k][1] = in->csignal[i][k][1] / prod[i][i];
			}
		}
		for (i = 0; i < NUM_OF_SIGNALS; i++) {
			for (j = 0; j < NUM_OF_SIGNALS; j++) {
				prod[i][j] = inner_Product_Old(in->csignal[i], in->csignal[j], in->psd, minfr,
						maxfr);
			}
		}
		for (i = 0; i < 2; i++) {
			for (k = 0; k < in->size; k++) {
				out->csignal[2 * i + 1][k][0] = in->csignal[2 * i + 1][k][0]
						- in->csignal[2 * i][k][0] * prod[2 * i][2 * i + 1] / prod[2 * i][2 * i];
				out->csignal[2 * i + 1][k][1] = in->csignal[2 * i + 1][k][1]
						- in->csignal[2 * i][k][1] * prod[2 * i][2 * i + 1] / prod[2 * i][2 * i];
			}
		}
	} else {
		// the orthonormalised signals overwrites the original ones
		// normalising
		for (i = 0; i < NUM_OF_SIGNALS; i++) {
			for (j = 0; j < NUM_OF_SIGNALS; j++) {
				prod[i][j] = sqrt(inner_Product_Old(in->csignal[i], in->csignal[j], in->psd, minfr,
						maxfr));
				if (i == j) {
					for (k = 0; k < in->size; k++) {
						in->csignal[i][k][0] /= prod[i][j];
						in->csignal[i][k][1] /= prod[i][j];
					}
				}
			}
		}
		// orthogonising
		for (i = 0; i < NUM_OF_SIGNALS; i++) {
			for (j = 0; j < NUM_OF_SIGNALS; j++) {
				prod[i][j] = inner_Product_Old(in->csignal[i], in->csignal[j], in->psd, minfr,
						maxfr);
			}
		}
		for (i = 0; i < 2; i++) {
			for (k = 0; k < in->size; k++) {
				in->csignal[2 * i + 1][k][0] -= in->csignal[2 * i][k][0] * prod[2 * i][2 * i + 1]
						/ prod[2 * i][2 * i];
				in->csignal[2 * i + 1][k][1] -= in->csignal[2 * i][k][1] * prod[2 * i][2 * i + 1]
						/ prod[2 * i][2 * i];
			}
		}
	}
}

double match_Typical(signalStruct *in, long minfr, long maxfr) {
	assert(in);
	assert(minfr>0);
	assert(minfr< maxfr);
	short i, j;
	for (i = 0; i < NUM_OF_SIGNALS; i++) {
		fftw_execute(in->plan[i]);
	}
	orthonormalise(in, minfr, maxfr, in);
	double prod[NUM_OF_SIGNALS][NUM_OF_SIGNALS];
	for (i = 0; i < NUM_OF_SIGNALS; i++) {
		for (j = 0; j < NUM_OF_SIGNALS; j++) {
			prod[i][j] = inner_Product_Old(in->csignal[i], in->csignal[j], in->psd, minfr, maxfr);
		}
	}
	long k;
	double M, Mmax = 0.;
	fftw_complex * new;
	new = fftw_malloc(in->size * sizeof(fftw_complex));
	fftw_plan iplan;
	for (i = 0; i < 2; i++) {
		iplan = fftw_plan_dft_c2r_1d(in->size, new, in->signal[i + 2], FFTW_ESTIMATE);
		cross_Product(in->csignal[H1P], in->csignal[i + 2], in->psd, minfr, maxfr, new);
		for (k = 0; k < in->size; k++) {
			if (k < minfr || k > maxfr) {
				new[k][0] = new[k][1] = 0.;
			}
		}
		fftw_execute(iplan);
		fftw_destroy_plan(iplan);
	}
	fftw_free(new);
	for (k = 0; k < in->size; k++) {// valamiért ki kell hagyni az első elemet
		M = sqrt(SQR(in->signal[H2P][k]) + SQR(in->signal[H2C][k]));
		Mmax = Mmax > M ? Mmax : M;
	}
	return Mmax / 2.; /// az IFFT f_min...f_max és -f_max...-f_min közötti tartományt is transzformálja
}

void calc_Overlap_Time(double *best, double *worst, signalStruct *in, long minfr, long maxfr) {
	assert(in);
	assert(minfr>0);
	assert(minfr< maxfr);
	short i, j;
	long k;
	for (i = 0; i < NUM_OF_SIGNALS; i++) {
		fftw_execute(in->plan[i]);
	}
	orthonormalise(in, minfr, maxfr, in);
	double prod[4][4];
	for (i = 0; i < NUM_OF_SIGNALS; i++) {
		for (j = 0; j < NUM_OF_SIGNALS; j++) {
			prod[i][j] = inner_Product_Old(in->csignal[i], in->csignal[j], in->psd, minfr, maxfr);
		}
	}
	///< \todo itt rossz
	fftw_complex *new = fftw_malloc(in->size * sizeof(fftw_complex));
	fftw_plan iplan;
	for (i = 0; i < NUM_OF_SIGNALS; i++) {
		cross_Product(in->csignal[i / 2], in->csignal[i % 2 + 2], in->psd, minfr, maxfr, new);
		for (k = 0; k < in->size; k++) {
			if (k < minfr || k > maxfr) {
				new[k][0] = new[k][1] = 0.;
			}
		}
		iplan = fftw_plan_dft_c2r_1d(in->size, new, in->signal[i], FFTW_ESTIMATE);
		fftw_execute(iplan);
		fftw_destroy_plan(iplan);
	}
	fftw_free(new);
	*best = match_Best(in, prod);
	*worst = match_Worst(in, prod);
}

double match_Best(signalStruct *in, double prod[NUM_OF_SIGNALS][NUM_OF_SIGNALS]) {
	assert(in);
	double A, B, C, M, maxfr_Match = 0.;
	long i;
	for (i = 0; i < in->size; i++) {
		A = SQR(in->signal[H1P][i]) + SQR(in->signal[H1C][i]);
		B = SQR(in->signal[H2P][i]) + SQR(in->signal[H2C][i]);
		C = in->signal[H1P][i] * in->signal[H2P][i] + in->signal[H1C][i] * in->signal[H2C][i];
		M = sqrt((A + B) / 2. + sqrt(SQR(A - B) / 4. + SQR(C)));
		maxfr_Match = maxfr_Match > M ? maxfr_Match : M;
	}
	return maxfr_Match / 2.; /// az IFFT f_min...f_max és -f_max...-f_min közötti tartományt is transzformálja
}

double match_Worst(signalStruct *in, double prod[NUM_OF_SIGNALS][NUM_OF_SIGNALS]) {
	assert(in);
	double A, B, C, M, maxfr_Match = 0.;
	long i;
	for (i = 0; i < in->size; i++) {
		A = SQR(in->signal[H1P][i]) + SQR(in->signal[H2P][i]);
		B = SQR(in->signal[H1C][i]) + SQR(in->signal[H2C][i]);
		C = in->signal[H1P][i] * in->signal[H2P][i] + in->signal[H1C][i] * in->signal[H2C][i];
		M = sqrt((A + B) / 2. - sqrt(SQR(A - B) / 4. + SQR(C)));
		maxfr_Match = maxfr_Match > M ? maxfr_Match : M;
	}
	return maxfr_Match / 2.; /// az IFFT f_min...f_max és -f_max...-f_min közötti tartományt is transzformálja
}
