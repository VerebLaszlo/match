/**
 * @file match.c
 * @author László Veréb
 * @date 2010.04.09.
 */

#include "match.h"

extern double pi;

int create_Signal_Struct(signalStruct *s, long size) {
	s->size = size;
	short i;
	for (i = 0; i < NUM_OF_SIGNALS; i++) {
		s->signal[i] = fftw_malloc(s->size * sizeof(double));
		s->csignal[i] = fftw_malloc(s->size * sizeof(fftw_complex));
		s->plan[i] = fftw_plan_dft_r2c_1d(s->size, s->signal[i], s->csignal[i],
				FFTW_ESTIMATE);
	}
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
		if (s->plan[i]) {
			fftw_destroy_plan(s->plan[i]);
		}
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
