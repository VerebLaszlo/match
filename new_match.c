/*
 * new_match.c
 *
 *  Created on: Sep 18, 2010
 *      Author: vereb
 */

#include <stdio.h>
#include <lal/LALDatatypes.h>		// LALStatus
#include <lal/GenerateInspiral.h>	// LALGenerateInspiral()
#include <lal/SimulateCoherentGW.h>	// CoherentGW
#include <lal/LIGOMetadataTables.h>		// SimInspiralTable
#include <lal/GeneratePPNInspiral.h>	// PPNParamStruc
#include <lal/LALSQTPNWaveformInterface.h>	// XLALSQTPNDestroyCoherentGW()
#include "match.h"

#define PREC "% -15.10lg "

double pi = LAL_PI;

typedef struct tagMatchStruct {
	double *signal[2];
	fftw_complex *csignal[2];
	fftw_plan plan[2];
	int length;
} matchStruct;

void mallocMatchStruct(matchStruct* m, int length) {
	m->length = length;
	m->signal[0] = fftw_malloc(m->length *sizeof(double));
	m->csignal[0] = fftw_malloc(m->length *sizeof(fftw_complex));
	m->plan[0] = fftw_plan_dft_r2c_1d(m->length, m->signal[0], m->csignal[0], FFTW_ESTIMATE);
	m->signal[1] = fftw_malloc(m->length *sizeof(double));
	m->csignal[1] = fftw_malloc(m->length *sizeof(fftw_complex));
	m->plan[1] = fftw_plan_dft_r2c_1d(m->length, m->signal[1], m->csignal[1], FFTW_ESTIMATE);
}

void freeMatchStruct(matchStruct *m) {
	fftw_free(m->signal[0]);
	fftw_free(m->csignal[0]);
	fftw_destroy_plan(m->plan[0]);
	fftw_free(m->signal[1]);
	fftw_free(m->csignal[1]);
	fftw_destroy_plan(m->plan[1]);
}

//****  PRÓBA  ****//
//****  PRÓBA  ****//

int main(int argc, char *argv[]) {
	static LALStatus status;
	CoherentGW waveform[2];
	SimInspiralTable injParams;
	PPNParamStruc ppnParams;
	matchStruct match;
	char *PNString[] = {"SpinQuadTaylortwoPointFivePNALL", "SpinQuadTaylortwoPointFivePNSS"};
	short num = 1;//12;	// number of the predefined configurations
	memset(&injParams, 0, sizeof(SimInspiralTable));
	memset(&ppnParams, 0, sizeof(PPNParamStruc));

	// kezdeti adatok beolvasása
	injParams.mass1 = atof(argv[1]);
	injParams.mass2 = atof(argv[2]);
	injParams.spin1x = atof(argv[3]);
	injParams.spin1y = atof(argv[4]);
	injParams.spin1z = atof(argv[5]);
	injParams.spin2x = atof(argv[6]);
	injParams.spin2y = atof(argv[7]);
	injParams.spin2z = atof(argv[8]);
	injParams.qmParameter1 = 1.;
	injParams.qmParameter2 = 1.;
	injParams.inclination = atof(argv[9]);
	injParams.f_lower = atof(argv[10]);
	injParams.distance = atof(argv[11]);
	ppnParams.deltaT = atof(argv[12]);
	injParams.polarization = 0;


	double freq_i, freq_f, freq_step, fr, dt;
	long index_i, index_f;
	freq_i = injParams.f_lower;
	dt = ppnParams.deltaT;
	// antenna függvények kiszámolása
	double t, p, s, fp, fc;
	t = p = s = 0.;
	fp = 0.5 * (1 + t * t) * cos(p) * cos(s) - t * sin(p) * sin(s);
	fc = 0.5 * (1 + t * t) * cos(p) * sin(s) + t * sin(p) * cos(s);
	double hp, hc, a1, a2, phi, shift;
	int i, length;
	short index, second, longer;
	for (index = 0; index < num; index++) {
		// hullámformák gyártása
		for (second = 0; second < 2; second++) {
			memset(&status, 0, sizeof(LALStatus));
			memset(&waveform[second], 0, sizeof(CoherentGW));
			LALSnprintf(injParams.waveform, LIGOMETA_WAVEFORM_MAX * sizeof(CHAR), PNString[second]);
			LALGenerateInspiral(&status, &waveform[second], &injParams, &ppnParams);
			if (status.statusCode) {
				fprintf(stderr, "LALSQTPNWaveformTest: error generating waveform\n");
				return status.statusCode;
			}
		}
		longer = waveform[0].f->data->length > waveform[1].f->data->length ? 0 : 1;
		length = waveform[longer].f->data->length;
		freq_f = waveform[longer].f->data->data[length - 1];
		mallocMatchStruct(&match, length);
		//****  PRÓBA  ****//
		FILE *file = fopen("wave.out", "w");
		//****  PRÓBA  ****//
		for (second = 0; second < 2; second++) {
			memset(match.signal[second], 0, match.length * sizeof(double));
			for (i = 0; i < waveform[second].f->data->length; i++) {
				a1 = waveform[second].a->data->data[2 * i];
				a2 = waveform[second].a->data->data[2 * i + 1];
				phi = waveform[second].phi->data->data[i] - waveform[second].phi->data->data[0];
				shift = waveform[second].shift->data->data[i];
				hp = a1 * cos(shift) * cos(phi) - a2 * sin(shift) * sin(phi);
				hc = a1 * sin(shift) * cos(phi) + a2 * cos(shift) * sin(phi);
				match.signal[second][i] = fp * hp + fc * hc;
				//****  PRÓBA  ****//
				fprintf(file, PREC PREC"\n", i*dt, match.signal[second][i]);fflush(file);
				//****  PRÓBA  ****//
			}
			XLALSQTPNDestroyCoherentGW(&waveform[second]);
			fftw_execute(match.plan[second]);
		}
		//****  PRÓBA  ****//
		fclose(file);
		//****  PRÓBA  ****//
		fr = 0.;
		freq_step = 1. / (dt * length);
		index_i = 0;
		while (fr < freq_i) {
			fr += freq_step;
			++index_i;
		}
		index_f = index_i;
		while (fr < freq_f) {
			fr += freq_step;
			index_f++;
		}
		freeMatchStruct(&match);
	}
	puts("Done.");
	return 0;
}
