/**
 * @file match.c
 * @author László Veréb
 * @date 2010.04.09.
 */

#include "match.h"

double pi;

int create_Signal_Struct(signalStruct *signal, long size) {
	signal->size = size;
	short i;
	for (i = 0; i < NUM_OF_SIGNALS; i++) {
		signal->signal[i] = fftw_malloc(signal->size * sizeof(double));
		memset(signal->signal[i], 0, signal->size * sizeof(double));
		signal->csignal[i] = fftw_malloc(signal->size * sizeof(fftw_complex));
		memset(signal->csignal[i], 0, signal->size * sizeof(fftw_complex));
		signal->plan[i] = fftw_plan_dft_r2c_1d(signal->size, signal->signal[i],
				signal->csignal[i], FFTW_ESTIMATE);
	}
	signal->psd = fftw_malloc(signal->size * sizeof(double));
	memset(signal->psd, 0, signal->size * sizeof(double));
	return MATCH_SUCCES;
}

void destroy_Signal_Struct(signalStruct *signal) {
	short i;
	for (i = 0; i < NUM_OF_SIGNALS; i++) {
		if (signal->signal[i]) {
			fftw_free(signal->signal[i]);
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

double inner_Product(fftw_complex left[], fftw_complex right[], double norm[],
		long minfr, long maxfr) {
	double scalar = 0.;
	long i;
	for (i = minfr; i < maxfr; i++) {
		scalar += ((left[i][0] * right[i][0] + left[i][1] * right[i][1])
				/ norm[i]);
	}
	return 4. * scalar;
}

void product(fftw_complex out[], fftw_complex left[], fftw_complex right[],
		double norm[], long minfr, long maxfr) {
	long i;
	for (i = minfr + 1; i < maxfr; i++) {
		out[i][0] = (left[i][0] * right[i][0] + left[i][1] * right[i][1])
				/ norm[i];
		out[i][1] = (left[i][1] * right[i][0] - left[i][0] * right[i][1])
				/ norm[i];
	}
}

void normalise(signalStruct *out, signalStruct *in, long minfr, long maxfr) {
	double prod;
	short i;
	long j;
	if (out) {
		// the normalised signals are stored in the out structure
		out->size = in->size;
		for (i = 0; i < NUM_OF_SIGNALS; i++) {
			prod = inner_Product(in->csignal[i], in->csignal[i], in->psd,
					minfr, maxfr);
			for (j = 0; j < in->size; j++) {
				out->csignal[i][j][0] = in->csignal[i][j][0] / prod;
				out->csignal[i][j][1] = in->csignal[i][j][1] / prod;
			}
		}
	} else {
		// the normalised signals overwrites the original ones
		for (i = 0; i < NUM_OF_SIGNALS; i++) {
			prod = sqrt(inner_Product(in->csignal[i], in->csignal[i], in->psd,
					minfr, maxfr));
			for (j = 0; j < in->size; j++) {
				in->csignal[i][j][0] /= prod;
				in->csignal[i][j][1] /= prod;
			}
		}
	}
}

void orthogonise(signalStruct *out, signalStruct *in, long minfr, long maxfr) {
	double prod[NUM_OF_SIGNALS][NUM_OF_SIGNALS];
	short i, j;
	long k;
	if (out) {
		// the orthonormalised signals are stored in the out structure
		out->size = in->size;
		double prod[NUM_OF_SIGNALS][NUM_OF_SIGNALS];
		for (i = 0; i < NUM_OF_SIGNALS; i++) {
			for (j = 0; j < NUM_OF_SIGNALS; j++) {
				prod[i][j] = inner_Product(in->csignal[i], in->csignal[j],
						in->psd, minfr, maxfr);
			}
		}
		for (i = 0; i < 2; i++) {
			for (k = 0; k < in->size; k++) {
				out->csignal[2 * i + 1][k][0] = in->csignal[2 * i + 1][k][0]
						- in->csignal[2 * i][k][0] * prod[2 * i][2 * i + 1]
								/ prod[2 * i][2 * i];
				out->csignal[2 * i + 1][k][1] = in->csignal[2 * i + 1][k][1]
						- in->csignal[2 * i][k][1] * prod[2 * i][2 * i + 1]
								/ prod[2 * i][2 * i];
			}
		}
	} else {
		// the orthonormalised signals overwrites the original ones
		for (i = 0; i < NUM_OF_SIGNALS; i++) {
			for (j = 0; j < NUM_OF_SIGNALS; j++) {
				prod[i][j] = inner_Product(in->csignal[i], in->csignal[j],
						in->psd, minfr, maxfr);
			}
		}
		for (i = 0; i < 2; i++) {
			for (k = 0; k < in->size; k++) {
				in->csignal[2 * i + 1][k][0] -= in->csignal[2 * i][k][0]
						* prod[2 * i][2 * i + 1] / prod[2 * i][2 * i];
				in->csignal[2 * i + 1][k][1] -= in->csignal[2 * i][k][1]
						* prod[2 * i][2 * i + 1] / prod[2 * i][2 * i];
			}
		}
	}
}

void orthonormalise(signalStruct *out, signalStruct *in, long minfr, long maxfr) {
	double prod[NUM_OF_SIGNALS][NUM_OF_SIGNALS];
	short i, j;
	long k;
	if (out) {
		// the orthonormalised signals are stored in the out structure
		out->size = in->size;
		for (i = 0; i < NUM_OF_SIGNALS; i++) {
			out->size = in->size;
			prod[i][i] = sqrt(inner_Product(in->csignal[i], in->csignal[i],
					in->psd, minfr, maxfr));
			for (k = 0; k < in->size; k++) {
				out->csignal[i][k][0] = in->csignal[i][k][0] / prod[i][i];
				out->csignal[i][k][1] = in->csignal[i][k][1] / prod[i][i];
			}
		}
		for (i = 0; i < NUM_OF_SIGNALS; i++) {
			for (j = 0; j < NUM_OF_SIGNALS; j++) {
				prod[i][j] = inner_Product(in->csignal[i], in->csignal[j],
						in->psd, minfr, maxfr);
			}
		}
		for (i = 0; i < 2; i++) {
			for (k = 0; k < in->size; k++) {
				out->csignal[2 * i + 1][k][0] = in->csignal[2 * i + 1][k][0]
						- in->csignal[2 * i][k][0] * prod[2 * i][2 * i + 1]
								/ prod[2 * i][2 * i];
				out->csignal[2 * i + 1][k][1] = in->csignal[2 * i + 1][k][1]
						- in->csignal[2 * i][k][1] * prod[2 * i][2 * i + 1]
								/ prod[2 * i][2 * i];
			}
		}
	} else {
		// the orthonormalised signals overwrites the original ones
		// normalising
		for (i = 0; i < NUM_OF_SIGNALS; i++) {
			for (j = 0; j < NUM_OF_SIGNALS; j++) {
				prod[i][j] = sqrt(inner_Product(in->csignal[i], in->csignal[j],
						in->psd, minfr, maxfr));
			}
			if (i == j) {
				for (k = 0; k < in->size; k++) {
					in->csignal[i][k][0] /= prod[i][j];
					in->csignal[i][k][1] /= prod[i][j];
				}
			}
		}
		// orthogonising
		for (i = 0; i < NUM_OF_SIGNALS; i++) {
			for (j = 0; j < NUM_OF_SIGNALS; j++) {
				prod[i][j] = inner_Product(in->csignal[i], in->csignal[j],
						in->psd, minfr, maxfr);
			}
		}
		for (i = 0; i < 2; i++) {
			for (k = 0; k < in->size; k++) {
				in->csignal[2 * i + 1][k][0] -= in->csignal[2 * i][k][0]
						* prod[2 * i][2 * i + 1] / prod[2 * i][2 * i];
				in->csignal[2 * i + 1][k][1] -= in->csignal[2 * i][k][1]
						* prod[2 * i][2 * i + 1] / prod[2 * i][2 * i];
			}
		}
	}
}

double match_Simple(signalStruct *in, long minfrfr, long maxfr) {
	short i;
	long j;
	for (i = 0; i < 2; i++) {
		for (j = 0; j < in->size; j++) {
			in->signal[2 * i][j] += in->signal[2 * i + 1][j];
		}
		fftw_execute(in->plan[2 * i]);
	}
	double prod12, prod11, prod22;
	prod12 = inner_Product(in->csignal[0], in->csignal[2], in->psd, minfrfr,
			maxfr);
	prod11 = inner_Product(in->csignal[0], in->csignal[0], in->psd, minfrfr,
			maxfr);
	prod22 = inner_Product(in->csignal[2], in->csignal[2], in->psd, minfrfr,
			maxfr);
	return prod12 / sqrt(prod11 * prod22);
}

double match_Typical(signalStruct *in, long minfr, long maxfr,
		maxMethods method) {
	short i, j;
	for (i = 0; i < NUM_OF_SIGNALS; i++) {
		fftw_execute(in->plan[i]);
	}
	orthonormalise(NULL, in, minfr, maxfr);
	double prod[NUM_OF_SIGNALS][NUM_OF_SIGNALS];
	for (i = 0; i < NUM_OF_SIGNALS; i++) {
		for (j = 0; j < NUM_OF_SIGNALS; j++) {
			prod[i][j] = inner_Product(in->csignal[i], in->csignal[j], in->psd,
					minfr, maxfr);
		}
	}
	if (method == TIME) {
		long k;
		double M, Mmax = 0.;
		fftw_complex * new;
		new = fftw_malloc(in->size * sizeof(fftw_complex));
		/*
		 for (k = 0; k < in->size; k++) {
		 new[k][0] = new[k][1] = 0.;
		 }
		 */
		fftw_plan iplan;
		for (i = 0; i < 2; i++) {
			//memset(new, 0, in->size *sizeof(fftw_complex));
			iplan = fftw_plan_dft_c2r_1d(in->size, new, in->signal[2 * i],
					FFTW_ESTIMATE);
			product(new, in->csignal[2 * i], in->csignal[2 * i + 1], in->psd,
					minfr, maxfr);
			fftw_execute(iplan);
			fftw_destroy_plan(iplan);
		}
		fftw_free(new);
		double norm = 4. * sqrt(sqrt(prod[H1P][H1P] * prod[H2P][H2P]
				* prod[H1C][H1C] * prod[H2C][H2C])); ///< \verify miért jó így????
		norm = 4. * sqrt(sqrt(prod[H1P][H1P] * prod[H2P][H2P] * prod[H1C][H1C]
				* prod[H2C][H2C])); ///< \verify miért jó így????
		for (k = 0; k < in->size; k++) {
			M = sqrt(SQR(in->signal[0][k]) / (prod[H1P][H1P] * prod[H1C][H1C])
					+ SQR(in->signal[2][k]) / (prod[H2P][H2P]
							* prod[H2C][H2C])) / 4.;
			Mmax = Mmax > M ? Mmax : M;
		}
		return Mmax;
	} else {
		return sqrt(SQR(prod[H1P][H2P]) / (prod[H1P][H1P] * prod[H2P][H2P])
				+ SQR(prod[H1P][H2C]) / (prod[H1P][H1P] * prod[H2C][H2C]));
	}
}

void calc_Overlap(double *best, double *worst, signalStruct *in, long minfr,
		long maxfr) {
	short i, j;
	for (i = 0; i < NUM_OF_SIGNALS; i++) {
		fftw_execute(in->plan[i]);
	}
	orthonormalise(NULL, in, minfr, maxfr);
	double prod[4][4];
	for (i = 0; i < NUM_OF_SIGNALS; i++) {
		for (j = 0; j < NUM_OF_SIGNALS; j++) {
			prod[i][j] = inner_Product(in->csignal[i], in->csignal[j], in->psd,
					minfr, maxfr);
		}
	}
	double A = SQR(prod[H1P][H2P]) / (prod[H1P][H1P] * prod[H2P][H2P])
			+ SQR(prod[H1P][H2C]) / (prod[H1P][H1P] * prod[H2C][H2C]);
	double B = SQR(prod[H1C][H2P]) / (prod[H1C][H1C] * prod[H2P][H2P])
			+ SQR(prod[H1C][H2C]) / (prod[H1C][H1C] * prod[H2C][H2C]);
	double C = (prod[H1P][H2P] * prod[H1C][H2P] / sqrt(prod[H1P][H1P]
			* prod[H1C][H1C] * SQR(prod[H2P][H2P])) + prod[H1P][H2C]
			* prod[H1C][H2C]) / sqrt(prod[H1P][H1P] * prod[H1C][H1C]
			* SQR(prod[H2C][H2C]));
	*best = sqrt((A + B) / 2. + sqrt(SQR(A - B) / 4. + SQR(C)));
	*worst = sqrt((A + B) / 2. - sqrt(SQR(A - B) / 4. + SQR(C)));
}

void calc_Overlap_Time(double *best, double *worst, signalStruct *in,
		long minfr, long maxfr) {
	short i, j;
	for (i = 0; i < NUM_OF_SIGNALS; i++) {
		fftw_execute(in->plan[i]);
	}
	orthonormalise(NULL, in, minfr, maxfr);
	double prod[4][4];
	for (i = 0; i < NUM_OF_SIGNALS; i++) {
		for (j = 0; j < NUM_OF_SIGNALS; j++) {
			prod[i][j] = inner_Product(in->csignal[i], in->csignal[j], in->psd,
					minfr, maxfr);
		}
	}
	///< \todo itt rossz
	fftw_complex *new = fftw_malloc(in->size * sizeof(fftw_complex));
	fftw_plan iplan;
	for (i = 0; i < NUM_OF_SIGNALS; i++) {
		product(new, in->csignal[i / 2], in->csignal[i % 2 + 2], in->psd,
				minfr, maxfr);
		iplan = fftw_plan_dft_c2r_1d(in->size, new, in->signal[i],
				FFTW_ESTIMATE);
		fftw_execute(iplan);
		fftw_destroy_plan(iplan);
	}
	fftw_free(new);
	*best = match_Best(in, prod);
	*worst = match_Worst(in, prod);
}

double match_Best(signalStruct *in, double prod[NUM_OF_SIGNALS][NUM_OF_SIGNALS]) {
	double A, B, C, M, maxfr_Match = 0.;
	long i;
	for (i = 0; i < in->size; i++) {
		//		printf("%ld %ld\n", i, in->size);fflush(stdout);
		A = SQR(in->signal[H1P][i]) / (prod[H1P][H1P] * prod[H2P][H2P])
				+ SQR(in->signal[H2P][i]) / (prod[H1P][H1P] * prod[H2C][H2C]);
		B = SQR(in->signal[H1C][i]) / (prod[H1C][H1C] * prod[H2P][H2P])
				+ SQR(in->signal[H2C][i]) / (prod[H1C][H1C] * prod[H2C][H2C]);
		C = in->signal[H1P][i] * in->signal[H2P][i] / sqrt(prod[H1P][H1P]
				* prod[H1C][H1C] * SQR(prod[H2P][H2P])) + in->signal[H1C][i]
				* in->signal[H2C][i] / sqrt(prod[H1P][H1P] * prod[H1C][H1C]
				* SQR(prod[H2C][H2C]));
		//		A /=128.;
		//		B /=128.;
		//		C/=128.;
		//		printf(PREC PREC PREC"\n", A+B, A-B, C);fflush(stdout);
		//		printf(PREC PREC PREC"\n", (A+B)/2., SQR(A-B)/4., SQR(C));fflush(stdout);
		//		printf(PREC PREC"\n", (A+B)/2., sqrt(SQR(A-B)/4.+ SQR(C)));fflush(stdout);
		M = sqrt((A + B) / 2. + sqrt(SQR(A - B) / 4. + SQR(C)));
		//		printf("0 = "PREC", "PREC" = "PREC", M = "PREC"\n", C, A, B, M);
		fflush(stdout);
		//		printf(PREC PREC PREC PREC"\n", A, B, C, M);fflush(stdout);
		maxfr_Match = maxfr_Match > M ? maxfr_Match : M;
	}
	return maxfr_Match;
}

double match_Worst(signalStruct *in,
		double prod[NUM_OF_SIGNALS][NUM_OF_SIGNALS]) {
	double A, B, C, M, maxfr_Match = 0.;
	long i;
	for (i = 0; i < in->size; i++) {
		A = SQR(in->signal[H1P][i]) / (prod[H1P][H1P] * prod[H2P][H2P])
				+ SQR(in->signal[H2P][i]) / (prod[H1P][H1P] * prod[H2C][H2C]);
		B = SQR(in->signal[H1C][i]) / (prod[H1C][H1C] * prod[H2P][H2P])
				+ SQR(in->signal[H2C][i]) / (prod[H1C][H1C] * prod[H2C][H2C]);
		C = in->signal[H1P][i] * in->signal[H2P][i] / sqrt(prod[H1P][H1P]
				* prod[H1C][H1C] * SQR(prod[H2P][H2P])) + in->signal[H1C][i]
				* in->signal[H2C][i] / sqrt(prod[H1P][H1P] * prod[H1C][H1C]
				* SQR(prod[H2C][H2C]));
		M = sqrt((A + B) / 2. - sqrt(SQR(A - B) / 4. + SQR(C)));
		//printf(PREC PREC"\n", (A + B) / 2., sqrt(SQR(A-B) / 4. + SQR(C)));
		fflush(stdout);
		/*		printf("M="PREC "M^2="PREC "="PREC"-("PREC")^2\n", sqrt((A + B) / 2.
		 - sqrt(SQR(A - B) / 4. + SQR(C))), (A + B) / 2. - sqrt(
		 SQR(A - B) / 4. + SQR(C)), (A + B) / 2., SQR(A - B) / 4.
		 + SQR(C));*/
		maxfr_Match = maxfr_Match < M ? maxfr_Match : M;
	}
	return maxfr_Match;
}

// Old Version Starts

void blackman(double array[], long length, double wn[]) {
	double w, swss = 0.;
	double n = length;
	long i;
	for (i = 0; i < length; i++) {
		w = 0.42 - 0.5 * cos(2. * pi * (double) i / n) + 0.08 * cos(4. * pi
				* (double) i / n);
		wn[i] = array[i] * w;
		swss += (w * w);
	}
	swss = sqrt(swss / n);
	for (i = 0; i < length; i++) {
		wn[i] /= swss;
	}
}

double* psd(double n[], long length, double dt, void(*winf)(double array[],
		long length, double wn[])) {
	winf(n, 0, n);
	double * wn = fftw_malloc(length * sizeof(double));
	fftw_complex *cn = fftw_malloc(length * sizeof(fftw_complex));
	fftw_plan pn = fftw_plan_dft_r2c_1d(length, wn, cn, FFTW_ESTIMATE);
	blackman(n, length, wn);
	fftw_execute(pn);
	long i;
	for (i = 0; i < length; i++) {
		wn[i] = (cn[i][0] * cn[i][0] + cn[i][1] * cn[i][1]) * 2. * dt / length;
	}
	fftw_destroy_plan(pn);
	fftw_free(cn);
	return wn;
}

void calc_Time_Corr(double h1[], double h2[], double dest[], long length) {
	long i, j;
	long x = length - 1;
	for (i = 0; i < 2 * length; i++) {
		dest[i] = 0.;
	}
	for (i = 0; i < length; i++) {
		for (j = 0; j < length - i; j++) {
			/*if (i < 5 && j < 5) {
			 printf("%ld %ld %ld %ld %ld\n", length,j, j + 1, x, i+length);
			 fflush(stdout);
			 printf("%lg\n",h1[j + 1]);
			 fflush(stdout);
			 printf("%lg\n",h2[j]);
			 fflush(stdout);
			 printf("%lg\n",dest[x]);
			 fflush(stdout);
			 printf("%lg\n",dest[i + length]);
			 fflush(stdout);
			 }*/
			dest[i + length] = (dest[x] += (h2[j] * h1[j + i]));
		}
		x--;
	}
}

void copy_Detector(detector_Struct source, detector_Struct *dest) {
	dest->length = source.length;
	dest->det = source.det;
	dest->t = source.t;
	dest->s = source.s;
	dest->n = source.n;
	dest->ct = source.ct;
	dest->cs = source.cs;
	dest->cn = source.cn;
	dest->pt = source.pt;
	dest->ps = source.ps;
	dest->pn = source.pn;
}

// Old Version Ends
