/**
 * @file lal_wrapper.c
 *
 * @date Apr 9, 2011
 * @author vereb
 */

#include <lal/LALSQTPNWaveformInterface.h>
#include <lal/LALNoiseModelsInspiral.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/GenerateInspiral.h>
#include <lal/LALInspiralBank.h>
#include <lal/LALDatatypes.h>
#include <lal/LALInspiral.h>
#include <lal/LALStdlib.h>
#include <lal/RealFFT.h>
#include "lal_wrapper.h"

int switchMode = LALSQTPN_PRECESSING;

/** LAL parameters.
 */
typedef struct LALParameters {
	LALStatus status; ///<a
	CoherentGW waveform[2]; ///<a
	SimInspiralTable injParams[2]; ///<a
	PPNParamStruc ppnParams[2]; ///<a
	RandomInspiralSignalIn randIn; ///<a
	short shorter; ///<a
	long min_Length; ///<a
	long max_Length; ///<a
	Approximant approx[2];
} LALParameters;

/**
 * Done
 * @param lalparams
 * @param parameters
 */
static void initLALParameters(LALParameters *lalparams, System_Parameters *parameters) {
	assert(lalparams);
	assert(parameters);
	memset(&lalparams->waveform, 0, 2 * sizeof(CoherentGW));
	memset(&lalparams->injParams, 0, 2 * sizeof(SimInspiralTable));
	memset(&lalparams->ppnParams, 0, 2 * sizeof(PPNParamStruc));
	memset(&lalparams->waveform, 0, 2 * sizeof(CoherentGW));
	memset(&lalparams->status, 0, sizeof(LALStatus));
	lalparams->ppnParams[0].deltaT = 1. / parameters->freq_Sampling;
	parameters->freq_Min = 40.;
	LALStatus status;
	memset(&status, 1, sizeof(LALStatus));
	for (short i = 0; i < 2; i++) {
		lalparams->injParams[i].mass1 = parameters->system[i].bh[0].m;
		lalparams->injParams[i].mass2 = parameters->system[i].bh[1].m;
		lalparams->injParams[i].spin1x = parameters->system[i].bh[0].chi[0];
		lalparams->injParams[i].spin1y = parameters->system[i].bh[0].chi[1];
		lalparams->injParams[i].spin1z = parameters->system[i].bh[0].chi[2];
		lalparams->injParams[i].spin2x = parameters->system[i].bh[1].chi[0];
		lalparams->injParams[i].spin2y = parameters->system[i].bh[1].chi[1];
		lalparams->injParams[i].spin2z = parameters->system[i].bh[1].chi[2];
		lalparams->injParams[i].inclination = parameters->system[i].incl;
		lalparams->injParams[i].f_lower = parameters->freq_Min;
		lalparams->injParams[i].distance = parameters->system[i].dist;
		lalparams->injParams[i].coa_phase = parameters->system[i].coaPhase = 0.;
		lalparams->injParams[i].f_lower = parameters->freq_Initial;
		lalparams->ppnParams[i].deltaT = 1. / parameters->freq_Sampling;
		lalparams->injParams[i].amp_order = parameters->amp_Code[i];
		snprintf(lalparams->injParams[i].waveform, LIGOMETA_WAVEFORM_MAX * sizeof(CHAR), "%s%s%s",
			parameters->approx[i], parameters->phase[i], parameters->spin[i]);
		if (strstr(parameters->approx[i], "SpinQuadTaylor")) {
			lalparams->approx[i] = SpinQuadTaylor;
		} else if (strstr(parameters->approx[i], "SpinTaylorFrameless")) {
			lalparams->approx[i] = SpinTaylorFrameless;
		} else if (strstr(parameters->approx[i], "SpinTaylor")) {
			lalparams->approx[i] = SpinTaylor;
		}
	}
}

/**
 * Done
 * @param lalparams
 * @param psd
 */
static void createPSD(LALParameters *lalparams, double *psd) {
	assert(lalparams);
	lalparams->randIn.powerSpectrumDensity.length = lalparams->max_Length;
	double df = 1. / lalparams->ppnParams[0].deltaT / lalparams->randIn.powerSpectrumDensity.length;
	lalparams->randIn.powerSpectrumDensity.data = (REAL8*) LALMalloc(sizeof(REAL8) * lalparams->randIn.powerSpectrumDensity.length);
	LALNoiseSpectralDensity(&lalparams->status, &lalparams->randIn.powerSpectrumDensity, &LALLIGOIPsd, df);
	for (unsigned long j = 0; j < lalparams->randIn.powerSpectrumDensity.length; j++) {
		psd[j] = lalparams->randIn.powerSpectrumDensity.data[j];
	}
}

char apr[2][FILENAME_MAX];

/** Sets the signal from CoherentGW's A1A2
 * @param i
 * @param sig
 * @param lal
 * @param F
 */
static void setSignal_From_A1A2(short i, signalStruct *sig, LALParameters *lal) {
	REAL8 a1, a2, phi, shift;
	for (long j = 0; j < sig->length[i]; j++) {
		a1 = lal->waveform[i].a->data->data[2 * j];
		a2 = lal->waveform[i].a->data->data[2 * j + 1];
		phi = lal->waveform[i].phi->data->data[j];
		shift = 0.0;
		sig->signal[2 * i][j] = a1 * cos(shift) * cos(phi) - a2 * sin(shift) * sin(phi);
		sig->signal[2 * i + 1][j] = a1 * sin(shift) * cos(phi) + a2 * cos(shift) * sin(phi);
	}
}

/** Sets the signal from CoherentGW's HPHC
 * @param i
 * @param sig
 * @param lal
 * @param F
 */
/*static void setSignal_From_HPHC(short i, signalStruct *sig, LALParameters *lal) {
 for (long j = 0; j < sig->length[i]; j++) {
 sig->signal[2 * i][j] = lal->waveform[i].h->data->data[2 * j];
 sig->signal[2 * i + 1][j] = lal->waveform[i].h->data->data[2 * j + 1];
 }
 }*/

/**
 *
 * @param i
 * @param signal
 * @param lal
 */
static void set_Signal_From_H(short i, signalStruct *signal, LALParameters *lal) {
	for (long j = 0; j < signal->length[i]; j++) {
		signal->signal[H1P + 2 * i][j] = lal->waveform[i].h->data->data[j];
		signal->signal[H1C + 2 * i][j] = lal->waveform[i].h->data->data[signal->length[i] + j];
	}
}

/**
 *
 * @param signal
 * @param lal
 */
static void create_Signal_Struct_From_LAL(signalStruct *signal, LALParameters *lal) {
	for (short i = 0; i < 2; i++) {
		switch (lal->approx[i]) {
		case SpinQuadTaylor:
		case SpinTaylorFrameless:
		case SpinTaylor:
			signal->length[i] = lal->ppnParams[i].length;
			break;
		default:
			break;
		}
	}
	signal->size = fmax(signal->length[0], signal->length[1]);
	create_Signal_Struct(signal, signal->size);
	for (short i = 0; i < 2; i++) {
		switch (lal->approx[i]) {
		case SpinQuadTaylor:
			setSignal_From_A1A2(i, signal, lal);
			break;
		case SpinTaylorFrameless:
			set_Signal_From_H(i, signal, lal);
			break;
		case SpinTaylor:
			setSignal_From_A1A2(i, signal, lal);
			break;
		default:
			break;
		}
	}
}

int generateWaveformPair(System_Parameters *parameters, signalStruct *signal) {
	static LALParameters lalparams;
	initLALParameters(&lalparams, parameters);
	for (short i = 0; i < 2; i++) {
		memset(&lalparams.status, 0, sizeof(LALStatus));
		LALGenerateInspiral(&lalparams.status, &lalparams.waveform[i], &lalparams.injParams[i],
			&lalparams.ppnParams[i]);
		if (lalparams.status.statusCode) {
			fprintf(stderr, "%d: LALSQTPNWaveformTest: error generating waveform\n", i);
			XLALSQTPNDestroyCoherentGW(&lalparams.waveform[0]);
			return NOT_FOUND;
		}
		parameters->system[i].coaPhase = lalparams.ppnParams[i].fStop;
		parameters->system[i].coaTime = lalparams.ppnParams[i].tc;
	}
	if (signal) {
		create_Signal_Struct_From_LAL(signal, &lalparams);
		createPSD(&lalparams, signal->psd);
		parameters->freq_Step = 1. / (lalparams.ppnParams[0].deltaT * signal->size);
	}
	XLALSQTPNDestroyCoherentGW(&lalparams.waveform[0]);
	XLALSQTPNDestroyCoherentGW(&lalparams.waveform[1]);
	XLALFree(lalparams.randIn.powerSpectrumDensity.data);
	LALCheckMemoryLeaks();
	return FOUND;
}
