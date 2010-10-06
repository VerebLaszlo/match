/**
 * @file match.c
 * @author László Veréb
 * @date 2010.04.09.
 */

#include "match.h"

double pi;

int create_Signal_Struct(signalStruct *s, long size) {
	s->size = size;
	short i;
	for (i = 0; i < NUM_OF_SIGNALS; i++) {
		s->signal[i] = fftw_malloc(s->size * sizeof(double));
		memset(s->signal[i], 0, s->size * sizeof(double));
		s->csignal[i] = fftw_malloc(s->size * sizeof(fftw_complex));
		memset(s->csignal[i], 0, s->size * sizeof(fftw_complex));
		s->cproduct[i] = fftw_malloc(s->size * sizeof(fftw_complex));
		memset(s->cproduct[i], 0, s->size * sizeof(fftw_complex));
		s->length[i] = s->size;
		s->plan[i] = fftw_plan_dft_r2c_1d(s->size, s->signal[i], s->csignal[i],
				FFTW_ESTIMATE);
		s->iplan[i] = fftw_plan_dft_c2r_1d(s->size, s->cproduct[i],
				s->signal[i], FFTW_ESTIMATE);
	}
	s->psd = fftw_malloc(s->size * sizeof(double));
	memset(s->psd, 0, s->size * sizeof(double));
	return MATCH_SUCCES;
}

void destroy_Signal_Struct(signalStruct *s) {
	short i;
	for (i = 0; i < NUM_OF_SIGNALS; i++) {
		if (s->signal[i]) {
			fftw_free(s->signal[i]);
		}
		if (s->csignal[i]) {
			fftw_free(s->csignal[i]);
		}
		if (s->cproduct[i]) {
			fftw_free(s->cproduct[i]);
		}
		if (s->plan[i]) {
			fftw_destroy_plan(s->plan[i]);
		}
		if (s->iplan[i]) {
			fftw_destroy_plan(s->iplan[i]);
		}
	}
	if (s->psd) {
		fftw_free(s->psd);
	}
}

double inner_Product(fftw_complex left[], fftw_complex right[], double norm[],
		long min, long max) {
	double scalar = 0.;
	long i;
	for (i = min; i < max; i++) {
		scalar += ((left[i][0] * right[i][0] + left[i][1] * right[i][1])
				/ norm[i]);
	}
	return 4. * scalar;
}

void normalise_Self(signalStruct *s, long min, long max) {
	double prod;
	short i;
	long j;
	for (i = 0; i < NUM_OF_SIGNALS; i++) {
		prod = inner_Product(s->csignal[i], s->csignal[i], s->psd, min, max);
		for (j = 0; j < s->length[i]; j++) {
			s->csignal[i][j][0] /= prod;
			s->csignal[i][j][1] /= prod;
		}
	}
}

void normalise(signalStruct *out, signalStruct *s, long min, long max) {
	/// \todo is it needed to copy everything from s to out?
	out->size = s->size;
	double prod;
	short i;
	long j;
	for (i = 0; i < NUM_OF_SIGNALS; i++) {
		out->length[i] = s->length[i];
		prod = inner_Product(s->csignal[i], s->csignal[i], s->psd, min, max);
		for (j = 0; j < s->length[i]; j++) {
			out->csignal[i][j][0] = s->csignal[i][j][0] / prod;
			out->csignal[i][j][1] = s->csignal[i][j][1] / prod;
		}
	}
}

void orthogonisate_Self(signalStruct *s, long min, long max) {
	double prod, temp;
	short i;
	long j;
	for (i = 0; i < 2; i++) {
		prod = inner_Product(s->csignal[2 * i], s->csignal[2 * i + 1], s->psd,
				min, max);
		temp = sqrt(1. - SQR(prod));
		for (j = 0; j < s->length[2 * i]; j++) {
			s->csignal[2 * i + 1][j][0] = (s->csignal[2 * i + 1][j][0]
					- s->csignal[2 * i][j][0] * prod) / temp;
			s->csignal[2 * i + 1][j][1] = (s->csignal[2 * i + 1][j][1]
					- s->csignal[2 * i][j][1] * prod) / temp;
		}
	}
}

void orthogonisate(signalStruct *out, signalStruct *s, long min, long max) {
	/// \todo is it needed to copy everything from s to out?
	out->size = s->size;
	double prod, temp;
	short i;
	long j;
	for (i = 0; i < 2; i++) {
		out->length[2 * i] = s->length[2 * i];
		out->length[2 * i + 1] = s->length[2 * i + 1];
		prod = inner_Product(s->csignal[2 * i], s->csignal[2 * i + 1], s->psd,
				min, max);
		temp = sqrt(1. - SQR(prod));
		for (j = 0; j < s->length[2 * i]; j++) {
			out->csignal[2 * i + 1][j][0] = (s->csignal[2 * i + 1][j][0]
					- s->csignal[2 * i][j][0] * prod) / temp;
			out->csignal[2 * i + 1][j][1] = (s->csignal[2 * i + 1][j][1]
					- s->csignal[2 * i][j][1] * prod) / temp;
		}
	}
}

void proba(signalStruct *s, long min, long max) {
	double prod[NUM_OF_SIGNALS][NUM_OF_SIGNALS];
	short i, j;
	long k;
	for (i = 0; i < NUM_OF_SIGNALS; i++) {
		fftw_execute(s->plan[i]);
	}
	for (i = 0; i < NUM_OF_SIGNALS; i++) {
		for (j = 0; j < NUM_OF_SIGNALS; j++) {
			prod[i][j] = inner_Product(s->csignal[i], s->csignal[j], s->psd,
					min, max);
		}
		if (i == j) {
			for (k = 0; k < s->length[i]; k++) {
				s->csignal[i][k][0] /= prod[i][j];
				s->csignal[i][k][1] /= prod[i][j];
			}
		}
	}
	for (i = 0; i < NUM_OF_SIGNALS; i++) {
		for (j = 0; j < NUM_OF_SIGNALS; j++) {
			prod[i][j] = inner_Product(s->csignal[i], s->csignal[j], s->psd,
					min, max);
		}
	}
	printf("CDCD"PREC"\n", inner_Product(s->csignal[1], s->csignal[0], s->psd,
			min, max));
	for (k = 0; k < s->size; k++) {
		s->csignal[1][k][0] -= s->csignal[0][k][0] * prod[0][1] / prod[0][0];
		s->csignal[1][k][1] -= s->csignal[0][k][1] * prod[0][1] / prod[0][0];
	}
	printf("CDCD"PREC"\n", inner_Product(s->csignal[1], s->csignal[0], s->psd,
			min, max));
	fflush(stdout);
	//	printf("R"PREC PREC"\n", inner_Product(s->csignal[1], s->csignal[0], s->psd, min, max), prod[0][1]);
}

void orthonormalise_Self(signalStruct *s, long min, long max) {
	double prod[NUM_OF_SIGNALS][NUM_OF_SIGNALS];
	short i, j;
	long k;
	for (i = 0; i < NUM_OF_SIGNALS; i++) {
		for (j = 0; j < NUM_OF_SIGNALS; j++) {
			prod[i][j] = inner_Product(s->csignal[i], s->csignal[j], s->psd,
					min, max);
		}
		if (i == j) {
			for (k = 0; k < s->length[i]; k++) {
				s->csignal[i][k][0] /= prod[i][j];
				s->csignal[i][k][1] /= prod[i][j];
			}
		}
	}
	for (i = 0; i < NUM_OF_SIGNALS; i++) {
		for (j = 0; j < NUM_OF_SIGNALS; j++) {
			prod[i][j] = inner_Product(s->csignal[i], s->csignal[j], s->psd,
					min, max);
		}
	}
	for (i = 0; i < 2; i++) {
		for (k = 0; k < s->size; k++) {
			s->csignal[2*i+1][k][0] -= s->csignal[2*i][k][0] * prod[2*i][2*i+1]
					/ prod[2*i][2*i];
			s->csignal[2*i+1][k][1] -= s->csignal[2*i][k][1] * prod[2*i][2*i+1]
					/ prod[2*i][2*i];
		}
	}
	/*
	 double prod1;
	 for (i = 0; i < 2; i++) {
	 prod1 = inner_Product(s->csignal[2 * i], s->csignal[2 * i + 1], s->psd,
	 min, max);
	 temp = sqrt(1. - SQR(prod1));
	 for (j = 0; j < s->length[2 * i]; j++) {
	 s->csignal[2 * i + 1][j][0] = (s->csignal[2 * i + 1][j][0]
	 - s->csignal[2 * i][j][0] * prod1) / temp;
	 s->csignal[2 * i + 1][j][1] = (s->csignal[2 * i + 1][j][1]
	 - s->csignal[2 * i][j][1] * prod1) / temp;
	 }
	 }*/
}

void orthonormalise(signalStruct *out, signalStruct *s, long min, long max) {
	/// \todo is it needed to copy everything from s to out?
	out->size = s->size;
	double prod, temp;
	short i;
	long j;
	for (i = 0; i < NUM_OF_SIGNALS; i++) {
		out->length[i] = s->length[i];
		prod = inner_Product(s->csignal[i], s->csignal[i], s->psd, min, max);
		for (j = 0; j < s->length[i]; j++) {
			out->csignal[i][j][0] = s->csignal[i][j][0] / prod;
			out->csignal[i][j][1] = s->csignal[i][j][1] / prod;
		}
	}
	for (i = 0; i < 2; i++) {
		out->length[2 * i] = s->length[2 * i];
		out->length[2 * i + 1] = s->length[2 * i + 1];
		prod = inner_Product(s->csignal[2 * i], s->csignal[2 * i + 1], s->psd,
				min, max);
		temp = sqrt(1. - SQR(prod));
		for (j = 0; j < s->length[2 * i]; j++) {
			out->csignal[2 * i + 1][j][0] = (s->csignal[2 * i + 1][j][0]
					- s->csignal[2 * i][j][0] * prod) / temp;
			out->csignal[2 * i + 1][j][1] = (s->csignal[2 * i + 1][j][1]
					- s->csignal[2 * i][j][1] * prod) / temp;
		}
	}
}

void calculate_Constants(double *A, double *B, double *C, signalStruct *s,
		long min, long max) {
	double pp = inner_Product(s->csignal[0], s->csignal[2], s->psd, min, max);
	double cc = inner_Product(s->csignal[1], s->csignal[3], s->psd, min, max);
	double pc = inner_Product(s->csignal[0], s->csignal[3], s->psd, min, max);
	double cp = inner_Product(s->csignal[1], s->csignal[2], s->psd, min, max);
	*A = SQR(pp) + SQR(pc);
	*B = SQR(cp) + SQR(cc);
	*C = pp * cp + pc * cc;
}

double match_simple(signalStruct *s, long min, long max) {
	short i;
	long j;
	for (i = 0; i < 2; i++) {
		for (j = 0; j < s->length[i]; j++) {
			s->signal[2 * i][j] += s->signal[2 * i + 1][j];
		}
		fftw_execute(s->plan[2 * i]);
	}
	double prod12, prod11, prod22;
	prod12 = inner_Product(s->csignal[0], s->csignal[2], s->psd, min, max);
	prod11 = inner_Product(s->csignal[0], s->csignal[0], s->psd, min, max);
	prod22 = inner_Product(s->csignal[2], s->csignal[2], s->psd, min, max);
	return prod12 / sqrt(prod11 * prod22);
}

double match_typical(signalStruct *s, long min, long max) {
	short i;
	long j;
	for (i = 0; i < NUM_OF_SIGNALS; i++) {
		fftw_execute(s->plan[i]);
	}
	orthonormalise_Self(s, min, max);
	// number is the signal, the index is the polarisation 0=+, 1=x
	double prod11, prod22[2], prod12[2];
	prod11 = inner_Product(s->csignal[0], s->csignal[0], s->psd, min, max);
	for (i = 0; i < 2; i++) {
		prod12[i] = inner_Product(s->csignal[0], s->csignal[i + 2], s->psd,
				min, max);
		prod22[i] = inner_Product(s->csignal[i + 2], s->csignal[i + 2], s->psd,
				min, max);
	}
	return sqrt(SQR(prod12[0]) / (prod11 * prod22[0]) + SQR(prod12[1])
			/ (prod11 * prod22[1]));
	//	return sqrt(SQR(prod12[0]) / (prod11*prod22[0])*SQR(prod12[1]*55) / (prod11*prod22[1]));
	// átírni innen
	/*	for (i = 0; i < 2; i++) {
	 for (j = 0; j < s->length[0]; j++) {
	 if (s->psd[j] != 0.) {
	 s->cproduct[i][j][0] = (s->csignal[i / 2][j][0] * s->csignal[i % 2
	 + 2][j][0] + s->csignal[i / 2][j][1]
	 * s->csignal[i % 2 + 2][j][1]) / s->psd[j];
	 s->cproduct[i][j][1] = (s->csignal[i / 2][j][1] * s->csignal[i % 2
	 + 2][j][0] - s->csignal[i / 2][j][0]
	 * s->csignal[i % 2 + 2][j][1]) / s->psd[j];
	 }
	 }
	 }
	 */// átírni eddig
	/*	double M, Mmax = 0.;
	 for (i = 0; i < 2; i++) {
	 fftw_execute(s->iplan[i]);
	 }
	 for (j = 0; j < s->length[0]; j++) {
	 printf("%ld"PREC "^2"PREC"," PREC"^2" PREC"\n", j, s->signal[0][j], SQR(s->signal[0][j]), s->signal[1][j], SQR(s->signal[1][j]));fflush(stdout);
	 M = sqrt(SQR(s->signal[0][j]) + SQR(s->signal[1][j]));
	 Mmax = Mmax > M ? Mmax : M;
	 }
	 return Mmax;*/
}

double match_typical_time(signalStruct *s, long min, long max) {
	short i;
	long j;
	for (i = 0; i < NUM_OF_SIGNALS; i++) {
		fftw_execute(s->plan[i]);
	}
	orthonormalise_Self(s, min, max);
	// number is the signal, the index is the polarisation 0=+, 1=x
	double prod11, prod22[2], prod12[2];
	prod11 = inner_Product(s->csignal[0], s->csignal[0], s->psd, min, max);
	for (i = 0; i < 2; i++) {
		prod12[i] = inner_Product(s->csignal[0], s->csignal[i + 2], s->psd,
				min, max);
		prod22[i] = inner_Product(s->csignal[i + 2], s->csignal[i + 2], s->psd,
				min, max);
	}
	return sqrt(SQR(prod12[0]) / (prod11 * prod22[0]) + SQR(prod12[1])
			/ (prod11 * prod22[1]));
	//	return sqrt(SQR(prod12[0]) / (prod11*prod22[0])*SQR(prod12[1]*55) / (prod11*prod22[1]));
	// átírni innen
	/*	for (i = 0; i < 2; i++) {
	 for (j = 0; j < s->length[0]; j++) {
	 if (s->psd[j] != 0.) {
	 s->cproduct[i][j][0] = (s->csignal[i / 2][j][0] * s->csignal[i % 2
	 + 2][j][0] + s->csignal[i / 2][j][1]
	 * s->csignal[i % 2 + 2][j][1]) / s->psd[j];
	 s->cproduct[i][j][1] = (s->csignal[i / 2][j][1] * s->csignal[i % 2
	 + 2][j][0] - s->csignal[i / 2][j][0]
	 * s->csignal[i % 2 + 2][j][1]) / s->psd[j];
	 }
	 }
	 }
	 */// átírni eddig
	/*	double M, Mmax = 0.;
	 for (i = 0; i < 2; i++) {
	 fftw_execute(s->iplan[i]);
	 }
	 for (j = 0; j < s->length[0]; j++) {
	 printf("%ld"PREC "^2"PREC"," PREC"^2" PREC"\n", j, s->signal[0][j], SQR(s->signal[0][j]), s->signal[1][j], SQR(s->signal[1][j]));fflush(stdout);
	 M = sqrt(SQR(s->signal[0][j]) + SQR(s->signal[1][j]));
	 Mmax = Mmax > M ? Mmax : M;
	 }
	 return Mmax;*/
}

double match_Best(signalStruct *s) {
	double A, B, C, M, max_Match = 0.;
	long i;
	for (i = 0; i < s->length[0]; i++) {
		A = SQR(s->signal[0][i]) + SQR(s->signal[2][i]);
		B = SQR(s->signal[1][i]) + SQR(s->signal[3][i]);
		C = s->signal[0][i] * s->signal[2][i] + s->signal[1][i]
				* s->signal[3][i];
		M = sqrt((A + B) / 2. + sqrt(SQR(A - B) / 4. + SQR(C)));
		/*		printf("XXXXXXXX "PREC PREC PREC PREC"\n", s->signal[0][i],
		 s->signal[1][i], s->signal[2][i], s->signal[3][i]);
		 fflush(stdout);
		 */
		break;
		max_Match = max_Match > M ? max_Match : M;
	}
	return max_Match;
}

double match_Worst(signalStruct *s) {
	double A, B, C, M, max_Match = 0.;
	long i;
	for (i = 0; i < s->length[0]; i++) {
		A = SQR(s->signal[0][i]) + SQR(s->signal[2][i]);
		B = SQR(s->signal[1][i]) + SQR(s->signal[3][i]);
		C = s->signal[0][i] * s->signal[2][i] + s->signal[1][i]
				* s->signal[3][i];
		M = sqrt((A + B) / 2. - sqrt(SQR(A - B) / 4. + SQR(C)));
		printf("M="PREC "M^2="PREC "="PREC"-("PREC")^2\n", sqrt((A + B) / 2.
				- sqrt(SQR(A - B) / 4. + SQR(C))), (A + B) / 2. - sqrt(
				SQR(A - B) / 4. + SQR(C)), (A + B) / 2., SQR(A - B) / 4.
				+ SQR(C));
		max_Match = max_Match > M ? max_Match : M;
	}
	return max_Match;
}

void calc_Overlap(double *best, double *worst, signalStruct *s, long min,
		long max) {
	short i;
	long j;
	for (i = 0; i < NUM_OF_SIGNALS; i++) {
		fftw_execute(s->plan[i]);
	}
	orthonormalise_Self(s, min, max);
	/*	for (i = 0; i < NUM_OF_SIGNALS; i++) {
	 for (j = 0; j < 2; j++) {
	 printf("TEST: "PREC PREC PREC PREC PREC PREC"\n", s->signal[i][j],
	 s->signal[i][j], s->csignal[i][j][0], s->csignal[i][j][1], s->cproduct[i][j][0], s->cproduct[i][j][1]);
	 }
	 }*/
	///< \todo itt rossz
	for (i = 0; i < NUM_OF_SIGNALS; i++) {
		for (j = 0; j < s->length[0]; j++) {
			if (s->psd[j] != 0.) {
				s->cproduct[i][j][0] = (s->csignal[i / 2][j][0] * s->csignal[i
						% 2 + 2][j][0] + s->csignal[i / 2][j][1] * s->csignal[i
						% 2 + 2][j][1]) / s->psd[j];
				s->cproduct[i][j][1] = (s->csignal[i / 2][j][1] * s->csignal[i
						% 2 + 2][j][0] - s->csignal[i / 2][j][0] * s->csignal[i
						% 2 + 2][j][1]) / s->psd[j];
			}
		}
	}
	for (i = 0; i < NUM_OF_SIGNALS; i++) {
		fftw_execute(s->iplan[i]);
	}
	*best = match_Best(s);
	*worst = match_Worst(s);
}

// Old Version Starts

void mallocMatchStruct(matchStruct* m, int length) {
	m->length = length;
	m->signal[0] = fftw_malloc(m->length * sizeof(double));
	m->csignal[0] = fftw_malloc(m->length * sizeof(fftw_complex));
	m->plan[0] = fftw_plan_dft_r2c_1d(m->length, m->signal[0], m->csignal[0],
			FFTW_ESTIMATE);
	m->signal[1] = fftw_malloc(m->length * sizeof(double));
	m->csignal[1] = fftw_malloc(m->length * sizeof(fftw_complex));
	m->plan[1] = fftw_plan_dft_r2c_1d(m->length, m->signal[1], m->csignal[1],
			FFTW_ESTIMATE);
}

void freeMatchStruct(matchStruct *m) {
	fftw_free(m->signal[0]);
	fftw_free(m->csignal[0]);
	fftw_destroy_plan(m->plan[0]);
	fftw_free(m->signal[1]);
	fftw_free(m->csignal[1]);
	fftw_destroy_plan(m->plan[1]);
}

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

void multi_Malloc(long len, detector_Struct * det) {
	det->length = len;
	det->t = fftw_malloc(det->length * sizeof(double));
	det->s = fftw_malloc(det->length * sizeof(double));
	det->n = fftw_malloc(det->length * sizeof(double));
	det->ct = fftw_malloc(det->length * sizeof(fftw_complex));
	det->cs = fftw_malloc(det->length * sizeof(fftw_complex));
	det->cn = fftw_malloc(det->length * sizeof(fftw_complex));
	det->pt = fftw_plan_dft_r2c_1d(det->length, det->t, det->ct, FFTW_ESTIMATE);
	det->ps = fftw_plan_dft_r2c_1d(det->length, det->s, det->cs, FFTW_ESTIMATE);
	det->pn = fftw_plan_dft_r2c_1d(det->length, det->n, det->cn, FFTW_ESTIMATE);
}

void multi_Free(detector_Struct *det) {
	fftw_free(det->t);
	fftw_free(det->s);
	fftw_free(det->n);
	fftw_free(det->ct);
	fftw_free(det->cs);
	fftw_free(det->cn);
	fftw_destroy_plan(det->pt);
	fftw_destroy_plan(det->ps);
	fftw_destroy_plan(det->pn);
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

void calculate_After(long index, detector_Struct *det[], long *len) {
	long p1 = *len - 1, p2 = *len - 2;
	long i;
	long data_length = (*det)[p1].length;
	detector_Struct *temp;
	if (!index) {
		long length = *len + 2;
		//	det = realloc(det, length);
		temp = malloc(length * sizeof(detector_Struct));
		for (i = 0; i < *len; i++) {
			copy_Detector((*det)[i], &temp[i]);
		}
		free(*det);
		*det = temp;
	} else {
		temp = *det;
	}
	// kivenni a függvény elé
	multi_Malloc(data_length, &temp[*len]); // deallocated line: 221
	fftw_plan ipn = fftw_plan_dft_c2r_1d(temp[*len].length, temp[*len].cn,
			temp[*len].n, FFTW_ESTIMATE); // deallocated line: 71
	for (i = 0; i < temp[*len].length; i++) {
		temp[*len].ct[i][0] = temp[p2].ct[i][0] * temp[p1].ct[i][0]
				+ temp[p2].ct[i][1] * temp[p1].ct[i][1];
		temp[*len].ct[i][1] = temp[p2].ct[i][0] * temp[p1].ct[i][1]
				- temp[p2].ct[i][1] * temp[p1].ct[i][0];
		temp[*len].cs[i][0] = temp[p2].cs[i][0] * temp[p1].cs[i][0]
				+ temp[p2].cs[i][1] * temp[p1].cs[i][1];
		temp[*len].cs[i][1] = temp[p2].cs[i][0] * temp[p1].cs[i][1]
				- temp[p2].cs[i][1] * temp[p1].cs[i][0];
		temp[*len].cn[i][0] = temp[p2].cn[i][0] * temp[p1].cn[i][0]
				+ temp[p2].cn[i][1] * temp[p1].cn[i][1];
		temp[*len].cn[i][1] = temp[p2].cn[i][0] * temp[p1].cn[i][1]
				- temp[p2].cn[i][1] * temp[p1].cn[i][0];
	}
	fftw_execute(ipn);
	fftw_destroy_plan(ipn);
	++*len;
	if (!index) {
		multi_Malloc(data_length * 2, &temp[*len]); // deallocated line: 221
		printf("Correlation: ");
		fflush(stdout);
		calc_Time_Corr(temp[p2].t, temp[p1].t, temp[*len].t, temp[*len].length
				/ 2);
		printf("h ");
		fflush(stdout);
		calc_Time_Corr(temp[p2].s, temp[p1].s, temp[*len].s, temp[*len].length
				/ 2);
		printf("s ");
		fflush(stdout);
		calc_Time_Corr(temp[p2].n, temp[p1].n, temp[*len].n, temp[*len].length
				/ 2);
		printf("n\n");
		fflush(stdout);
		fftw_execute(temp[*len].pt);
		fftw_execute(temp[*len].ps);
		++*len;
	}
}

// Old Version Ends
