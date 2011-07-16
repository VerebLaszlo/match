/**
 * @file lal_wrapper.c
 *
 * @date Apr 9, 2011
 * @author vereb
 */

#include "lal_wrapper.h"

void createPSD(LALParameters *lalparams, double *psd) {
	assert(lalparams);
	lalparams->randIn.psd.length = lalparams->max_Length;
	double df = 1. / lalparams->ppnParams->deltaT / lalparams->randIn.psd.length;
	lalparams->randIn.psd.data = (REAL8*)LALMalloc(sizeof(REAL8) * lalparams->randIn.psd.length);
	LALNoiseSpectralDensity(&lalparams->status, &lalparams->randIn.psd, &LALLIGOIPsd, df);
	for (unsigned long j = 0; j < lalparams->randIn.psd.length; j++) {
		psd[j] = lalparams->randIn.psd.data[j];
	}
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
