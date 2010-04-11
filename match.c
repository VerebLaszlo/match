/**
 * @file match.c
 * @author László Veréb
 * @date 2010.04.09.
 */

#include <time.h>
#include "util.h"
#include "match.h"
#include "generator.h"
#include "detector.h"
#include "ioHandling.h"

extern double pi;

double scalar_freq(fftw_complex left[], fftw_complex right[], double norm[], long min, long max) {
	double scalar = 0.;
	long i;
	for (i = min; i < max; i++) {
		scalar += ((left[i][0] * right[i][0] + left[i][1] * right[i][1]) / norm[i]);
	}
	return 4. * scalar;
}

void blackman(double array[], long length, double wn[]) {
	double w, swss = 0.;
	double n = length;
	long i;
	for (i = 0; i < length; i++) {
		w = 0.42 - 0.5 * cos(2. * pi * (double) i / n) + 0.08 * cos(4. * pi * (double) i / n);
		wn[i] = array[i] * w;
		swss += (w * w);
	}
	swss = sqrt(swss / n);
	for (i = 0; i < length; i++) {
		wn[i] /= swss;
	}
}

double* psd(double n[], long length, double dt, void(*winf)(double array[], long length, double wn[])) {
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
/**
 *	Időbeli korrelációt végrehajtó függvény.
 */
void calc_Time_Corr(double h1[], double h2[], double dest[], long length) {
	long i, j;
	long x=length-1;
	for (i = 0; i < 2*length;i++) {
		dest[i] = 0.;
	}
	for (i = 0; i < length; i++) {
		for (j = 0; j < length - i; j++) {
			dest[i+length] = (dest[x] += (h2[j] * h1[j + i]));
		}
		x--;
	}
}

void multi_Malloc(size_t len, detector_Struct * det) {
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

void multi_Free(detector_Struct det) {
	fftw_free(det.t);
	fftw_free(det.s);
	fftw_free(det.n);
	fftw_free(det.ct);
	fftw_free(det.cs);
	fftw_free(det.cn);
	fftw_destroy_plan(det.pt);
	fftw_destroy_plan(det.ps);
	fftw_destroy_plan(det.pn);
}

void calculate_After(size_t index, detector_Struct det[], size_t *len) {
	size_t p1 = *len - 1, p2 = *len - 2;
	size_t i, j;
	size_t length = (++*len) + 1;
	det = realloc(det, length);
	// kivenni a függvény elé
	multi_Malloc(det[p1].length, &det[*len]);										// deallocated line: 221
	fftw_plan ipn = fftw_plan_dft_c2r_1d(det[*len].length, det[*len].cn, det[*len].n, FFTW_ESTIMATE);	// deallocated line: 71
	for(j = 0; j < det[*len].length; j++) {
		det[*len].ct[j][0] = det[p2].ct[j][0] * det[p1].ct[j][0] + det[p2].ct[j][1] * det[p1].ct[j][1];
		det[*len].ct[j][1] = det[p2].ct[j][0] * det[p1].ct[j][1] - det[p2].ct[j][1] * det[p1].ct[j][0];
		det[*len].cs[j][0] = det[p2].cs[j][0] * det[p1].cs[j][0] + det[p2].cs[j][1] * det[p1].cs[j][1];
		det[*len].cs[j][1] = det[p2].cs[j][0] * det[p1].cs[j][1] - det[p2].cs[j][1] * det[p1].cs[j][0];
		det[*len].cn[j][0] = det[p2].cn[j][0] * det[p1].cn[j][0] + det[p2].cn[j][1] * det[p1].cn[j][1];
		det[*len].cn[j][1] = det[p2].cn[j][0] * det[p1].cn[j][1] - det[p2].cn[j][1] * det[p1].cn[j][0];
	}
	fftw_execute(ipn);
	fftw_destroy_plan(ipn);
	if(!index) {
		++*len;
		multi_Malloc(det[p1].length * 2, &det[*len]);								// deallocated line: 221
		printf("Correlation: ");
		fflush(stdout);
		calc_Time_Corr(det[p2].t, det[p1].t, det[*len].t, det[*len].length / 2);
		printf("h ");
		fflush(stdout);
		calc_Time_Corr(det[p2].s, det[p1].s, det[*len].s, det[*len].length / 2);
		printf("s ");
		fflush(stdout);
		calc_Time_Corr(det[p2].n, det[p1].n, det[*len].n, det[*len].length / 2);
		printf("n\n");
		fflush(stdout);
		fftw_execute(det[*len].pt);
		fftw_execute(det[*len].ps);
	}
}

