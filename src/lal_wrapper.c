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
	REAL8Vector psd;
	short shorter; ///<a
	long min_Length; ///<a
	long max_Length; ///<a
	Approximant approx[2];
} LALParameters;

static INT4 convertAmplitudeOrderFromString(const char *amplitudeOrder) {
	INT4 amplitudeCode;
	if (!strstr(amplitudeOrder, "00PN")) {
		amplitudeCode = 0;
	} else if (!strstr(amplitudeOrder, "05PN")) {
		amplitudeCode = 1;
	} else if (!strstr(amplitudeOrder, "10PN")) {
		amplitudeCode = 2;
	} else {
		amplitudeCode = -1;
	}
	return amplitudeCode;
}

/**
 * Done
 * @param lalparams
 * @param parameters
 */
static void initLALParameters(LALParameters *lalparams, SystemParameter *parameters) {
	assert(lalparams);
	assert(parameters);
	memset(&lalparams->waveform, 0, 2 * sizeof(CoherentGW));
	memset(&lalparams->injParams, 0, 2 * sizeof(SimInspiralTable));
	memset(&lalparams->ppnParams, 0, 2 * sizeof(PPNParamStruc));
	memset(&lalparams->waveform, 0, 2 * sizeof(CoherentGW));
	memset(&lalparams->status, 0, sizeof(LALStatus));
	lalparams->ppnParams[0].deltaT = 1. / parameters->samplingFrequency;
	LALStatus status;
	memset(&status, 1, sizeof(LALStatus));
	for (short i = 0; i < 2; i++) {
		lalparams->injParams[i].mass1 = parameters->system[i].mass.mass[0];
		lalparams->injParams[i].mass2 = parameters->system[i].mass.mass[0];
		lalparams->injParams[i].spin1x = parameters->system[i].spin[0].component[FIXED][X];
		lalparams->injParams[i].spin1y = parameters->system[i].spin[0].component[FIXED][Y];
		lalparams->injParams[i].spin1z = parameters->system[i].spin[0].component[FIXED][Z];
		lalparams->injParams[i].spin2x = parameters->system[i].spin[1].component[FIXED][X];
		lalparams->injParams[i].spin2y = parameters->system[i].spin[1].component[FIXED][Y];
		lalparams->injParams[i].spin2z = parameters->system[i].spin[1].component[FIXED][Z];
		lalparams->injParams[i].inclination = parameters->system[i].inclination;
		lalparams->injParams[i].f_lower = parameters->initialFrequency;
		lalparams->injParams[i].distance = parameters->system[i].distance;
		//lalparams->injParams[i].coa_phase = parameters->system[i].coaPhase = 0.;
		lalparams->ppnParams[i].deltaT = 1. / parameters->samplingFrequency;
		lalparams->injParams[i].amp_order = convertAmplitudeOrderFromString(
			parameters->amplitude[i]);
		snprintf(lalparams->injParams[i].waveform, LIGOMETA_WAVEFORM_MAX * sizeof(CHAR), "%s%s%s",
			parameters->approximant[i], parameters->phase[i], parameters->spin[i]);
		if (strstr(parameters->approximant[i], "SpinQuadTaylor")) {
			lalparams->approx[i] = SpinQuadTaylor;
		} else if (strstr(parameters->approximant[i], "SpinTaylorFrameless")) {
			lalparams->approx[i] = SpinTaylorFrameless;
		} else if (strstr(parameters->approximant[i], "SpinTaylor")) {
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
	lalparams->psd.length = lalparams->max_Length;
	double df = 1. / lalparams->ppnParams[0].deltaT / lalparams->psd.length;
	lalparams->psd.data = (REAL8*) XLALCalloc(sizeof(REAL8), lalparams->psd.length);
	LALNoiseSpectralDensity(&lalparams->status, &lalparams->psd, &LALLIGOIPsd, df);
	for (unsigned long j = 0; j < lalparams->psd.length; j++) {
		psd[j] = lalparams->psd.data[j];
	}
}

char apr[2][FILENAME_MAX];

/** Sets the signal from CoherentGW's A1A2
 * @param i
 * @param sig
 * @param lal
 * @param F
 */
static void setSignal_From_A1A2(short i, SignalStruct *sig, LALParameters *lal) {
	REAL8 a1, a2, phi, shift;
	for (ulong j = 0; j < sig->length[i]; j++) {
		a1 = lal->waveform[i].a->data->data[2 * j];
		a2 = lal->waveform[i].a->data->data[2 * j + 1];
		phi = lal->waveform[i].phi->data->data[j];
		shift = 0.0;
		sig->componentsInTime[2 * i][j] = a1 * cos(shift) * cos(phi) - a2 * sin(shift) * sin(phi);
		sig->componentsInTime[2 * i + 1][j] = a1 * sin(shift) * cos(phi)
			+ a2 * cos(shift) * sin(phi);
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
static void set_Signal_From_H(short i, SignalStruct *signal, LALParameters *lal) {
	for (ulong j = 0; j < signal->length[i]; j++) {
		signal->componentsInTime[H1P + 2 * i][j] = lal->waveform[i].h->data->data[j];
		signal->componentsInTime[H1C + 2 * i][j] = lal->waveform[i].h->data->data[signal->length[i]
			+ j];
	}
}

/**
 *
 * @param signal
 * @param lal
 */
static void create_Signal_Struct_From_LAL(SignalStruct *signal, LALParameters *lal) {
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

int generateWaveformPair(SystemParameter *parameters, SignalStruct *signal) {
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
		parameters->coalescencePhase[i] = lalparams.ppnParams[i].fStop;
		parameters->coalescenceTime[i] = lalparams.ppnParams[i].tc;
	}
	if (signal) {
		create_Signal_Struct_From_LAL(signal, &lalparams);
		createPSD(&lalparams, signal->powerSpectrumDensity);
		parameters->samplingFrequency = 1. / (lalparams.ppnParams[0].deltaT * signal->size);
	}
	XLALSQTPNDestroyCoherentGW(&lalparams.waveform[0]);
	XLALSQTPNDestroyCoherentGW(&lalparams.waveform[1]);
	XLALFree(lalparams.psd.data);
	LALCheckMemoryLeaks();
	return FOUND;
}
