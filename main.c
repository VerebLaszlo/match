#include "util.h"
#include <lal/LALInspiralBank.h>
#include <lal/LALNoiseModelsInspiral.h>
#include <lal/RealFFT.h>
#include <lal/GenerateInspiral.h>
#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/LALSQTPNWaveformInterface.h>
#include "match.h"

double pi = M_PI;

long max_length = 60000;

int main(int argc, char *argv[]) {
	static LALStatus status;
	CoherentGW thewaveform1;
	CoherentGW thewaveform2;
	SimInspiralTable injParams;
	PPNParamStruc ppnParams;
	const char *filename;
	FILE *outputfile;
	INT4 i, length;
	char PNString[50];
	char SpinInter[50];

	if (argc != 18) {
			printf(
					"                         1  2  3   4   5   6   7   8   9    10		 11		 12		  13 14      	15			16			17\n");
			printf(
					"Correct parameter order: m1 m2 S1x S1y S1z S2x S2y S2z incl f_lower f_final distance dt PNorder1	PNorder2	SpinInter1	SpinInter2\n");
			return (1);
	}
	
	sprintf(PNString, "SpinQuadTaylor%s%s", argv[14], argv[16]);
	memset(&status, 0, sizeof(LALStatus));
	memset(&thewaveform1, 0, sizeof(CoherentGW));
	memset(&thewaveform2, 0, sizeof(CoherentGW));
	memset(&injParams, 0, sizeof(SimInspiralTable));
	memset(&ppnParams, 0, sizeof(PPNParamStruc));

	//	setting the parameters
	double noise[max_length];
	double time[max_length];
	static RandomInspiralSignalIn randIn;
	short j;
	size_t noise_length, max_Length;
	size_t diff1, diff2, diffn;
	FILE *file = fopen("noise.wave", "r");
	for (noise_length = 0; !feof(file) && noise_length < max_length; noise_length++) {
		i = fscanf(file, "%lg %lg\n", &time[noise_length], &noise[noise_length]);
	}
	fclose(file);
	double dt = time[1] - time[0];
	double t, p, s;
	t = p = s = 0.;
	injParams.mass1 = atof(argv[1]);
	injParams.mass2 = atof(argv[2]);
	injParams.spin1x = atof(argv[3]);
	injParams.spin1y = atof(argv[4]);
	injParams.spin1z = atof(argv[5]);
	injParams.spin2x = atof(argv[6]);
	injParams.spin2y = atof(argv[7]);
	injParams.spin2z = atof(argv[8]);
	injParams.qmParameter1 = 1.;//atof(argv[9]);
	injParams.qmParameter2 = 1.;//atof(argv[10]);
	injParams.inclination = atof(argv[9]);
	double freq_Min = injParams.f_lower = atof(argv[10]);
	double freq_Max = injParams.f_final = atof(argv[11]);
	injParams.distance = atof(argv[12]);
	ppnParams.deltaT = dt;// atof(argv[13]);
	injParams.polarization = 0;
	LALSnprintf(injParams.waveform, LIGOMETA_WAVEFORM_MAX * sizeof(CHAR), PNString);

	//interface(&mystatus, &thewaveform, &injParams, &ppnParams);
	LALGenerateInspiral(&status, &thewaveform1, &injParams, &ppnParams);
	if (status.statusCode) {
		fprintf( stderr, "LALSQTPNWaveformTest: error generating waveform\n" );
		return status.statusCode;
	}
	sprintf(PNString, "SpinQuadTaylor%s%s", argv[15], argv[17]);
	LALSnprintf(injParams.waveform, LIGOMETA_WAVEFORM_MAX * sizeof(CHAR), PNString);
	//interface(&mystatus, &thewaveform, &injParams, &ppnParams);
	LALGenerateInspiral(&status, &thewaveform2, &injParams, &ppnParams);
	if (status.statusCode) {
		fprintf( stderr, "LALSQTPNWaveformTest: error generating waveform\n" );
		return status.statusCode;
	}

	max_Length = thewaveform1.f->data->length > thewaveform2.f->data->length ? thewaveform1.f->data->length : thewaveform2.f->data->length;
	randIn.psd.length = max_Length;
	double df = 1. / dt / randIn.psd.length;
	randIn.psd.data = (REAL8*) LALMalloc(sizeof(REAL8)*randIn.psd.length);
	detector_Struct det_old, det_new;
	multi_Malloc(max_Length, &det_old);
	multi_Malloc(max_Length, &det_new);
	double hp, hc, a1, a2, phi, shift;
	double fp = 0.5*(1 + t*t)*cos(p)*cos(s) - t*sin(p)*sin(s);
    double fc = 0.5*(1 + t*t)*cos(p)*sin(s) + t*sin(p)*cos(s);
	double *norm;
	double *wn;
	short k;
	for (j = 0; j < 2; j++) {
		k = (j + 1) % 2;
		diff1 = (max_Length - thewaveform1.f->data->length);
		diff2 = (max_Length - thewaveform2.f->data->length);
		for (i = 0; i < max_Length; i++) {
			det_old.n[i] = det_new.n[i] = noise[noise_length - i];
		}
		for (i = 0; i < j * diff1; i++) {
			det_old.t[i] = det_new.t[i] = 0.;
		}
		for (; i < max_Length - k * diff1; i++) {
			a1  = thewaveform1.a->data->data[2*i];
			a2  = thewaveform1.a->data->data[2*i+1];
			phi     = thewaveform1.phi->data->data[i] - thewaveform1.phi->data->data[0];
			shift   = thewaveform1.shift->data->data[i];
			hp = a1*cos(shift)*cos(phi) - a2*sin(shift)*sin(phi);
			hc = a1*sin(shift)*cos(phi) + a2*cos(shift)*sin(phi);
			det_new.t[i] = fp * hp + fc * hc;
			det_old.t[i] = det_new.t[i];
		}
		for (; i < max_Length; i++) {
			det_old.t[i] = det_new.t[i] = 0.;
		}
		for (i = 0; i < j * diff2; i++) {
			det_old.s[i] = det_new.s[i] = 0.;
		}
		for (; i < max_Length - k * diff2; i++) {
			a1  = thewaveform1.a->data->data[2*i];
			a2  = thewaveform1.a->data->data[2*i+1];
			phi     = thewaveform1.phi->data->data[i] - thewaveform1.phi->data->data[0];
			shift   = thewaveform1.shift->data->data[i];
			hp = a1*cos(shift)*cos(phi) - a2*sin(shift)*sin(phi);
			hc = a1*sin(shift)*cos(phi) + a2*cos(shift)*sin(phi);
			det_new.s[i] = fp * hp + fc * hc;
			det_old.s[i] = det_new.s[i];
		}
		for (; i < max_Length; i++) {
			det_old.s[i] = det_new.s[i] = 0.;
		}
		fftw_execute(det_old.pt);
		fftw_execute(det_old.ps);
		fftw_execute(det_new.pt);
		fftw_execute(det_new.ps);
		norm = psd(det_old.n, det_old.length, dt, blackman);
		LALNoiseSpectralDensity (&status, &randIn.psd, &LALLIGOIPsd, df);
		wn = fftw_malloc(max_Length * sizeof(double));
		for (i = 0; i < max_Length; i++) {
			wn[i] = randIn.psd.data[i];
		}
		double freq_Step, fr = 0.;
		freq_Step = 1. / (dt * max_Length);
		long minfr = 0, maxfr = 0;
		while (fr < freq_Min) {
			fr += freq_Step;
			maxfr = ++minfr;
		}
		while (fr < freq_Max) {
			fr += freq_Step;
			maxfr++;
		}
		double ts_old = scalar_freq(det_old.ct, det_old.cs, norm, minfr, maxfr);
		double tt_old = scalar_freq(det_old.ct, det_old.ct, norm, minfr, maxfr);
		double ss_old = scalar_freq(det_old.cs, det_old.cs, norm, minfr, maxfr);
		double ts_new = scalar_freq(det_new.ct, det_new.cs, wn, minfr, maxfr);
		double tt_new = scalar_freq(det_new.ct, det_new.ct, wn, minfr, maxfr);
		double ss_new = scalar_freq(det_new.cs, det_new.cs, wn, minfr, maxfr);
		printf("old: %lg, new: %lg\n", ts_old / sqrt(tt_old * ss_old), ts_new / sqrt(tt_new * ss_new));
		fftw_free(wn);
		fftw_free(norm);
	}

	XLALSQTPNDestroyCoherentGW(&thewaveform1);
	XLALSQTPNDestroyCoherentGW(&thewaveform2);
	multi_Free(&det_old);
	multi_Free(&det_new);

	return 0;
}
