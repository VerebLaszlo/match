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

#define SQR(a) ((a)*(a))

long max_length = 60000;

int main(int argc, char *argv[]) {
	static LALStatus status;
	CoherentGW thewaveform1;
	CoherentGW thewaveform2;
	SimInspiralTable injParams;
	PPNParamStruc ppnParams;
	INT4 i;
	char PNString[50];

	if (argc != 15) {
		printf(
				"                         1  2  3   4   5   6   7   8   9    10		 11		  12 13      14\n");
		printf(
				"Correct parameter order: m1 m2 S1x S1y S1z S2x S2y S2z incl f_lower distance dt PNorder SpinInter1\n");
		return (1);
	}

	memset(&status, 0, sizeof(LALStatus));
	memset(&thewaveform1, 0, sizeof(CoherentGW));
	memset(&thewaveform2, 0, sizeof(CoherentGW));
	memset(&injParams, 0, sizeof(SimInspiralTable));
	memset(&ppnParams, 0, sizeof(PPNParamStruc));

	//	semax_tting the parameters
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
	double incl = injParams.inclination = atof(argv[9]);
	double freq_Min = injParams.f_lower = atof(argv[10]);
	injParams.distance = atof(argv[11]);
	double dt = ppnParams.deltaT = atof(argv[12]);
	injParams.polarization = 0;
	double freq_Max = 600.;
	// változó paraméterek
	double chi[2] = { sqrt(SQR(injParams.spin1x) + SQR(injParams.spin1y) + SQR(
			injParams.spin1z)), sqrt(SQR(injParams.spin2x) + SQR(
			injParams.spin2y) + SQR(injParams.spin2z)) };
	double where[12][2][3] = { { { injParams.spin1x / chi[0], injParams.spin1y
			/ chi[0], injParams.spin1z / chi[0] }, { injParams.spin2x / chi[1],
			injParams.spin2y / chi[1], injParams.spin2z / chi[1] } },
			{ { sin(incl), 0., cos(incl) }, { sin(incl), 0., cos(incl) } }, //1  L_N ~  S_1 ~  S_2
			{ { -sin(incl), 0., -cos(incl) }, { -sin(incl), 0., -cos(incl) } }, //2 -L_N ~  S_1 ~  S_2
			{ { sin(incl), 0., cos(incl) }, { -sin(incl), 0., -cos(incl) } }, //3  L_N ~  S_1 ~ -S_2
			{ { -sin(incl), 0., -cos(incl) }, { sin(incl), 0., cos(incl) } }, //4  L_N ~ -S_1 ~  S_2
			{ { cos(incl), 0., -sin(incl) }, { cos(incl), 0., -sin(incl) } }, //5         S_1 ~  S_2
			{ { cos(incl), 0., -sin(incl) }, { -cos(incl), 0., sin(incl) } }, //6         S_1 ~ -S_2
			{ { 0., 1., 0. }, { 0., 1., 0. } }, // pályasíkos párhuzamos, egyező irány, 	VERY GOOD
			{ { 0., 1., 0. }, { 0., -1., 0. } }, // pályasíkos párhuzamos, különböző irány
			{ { cos(incl), 0., sin(incl) }, { 0., 1., 0. } }, // pályasíkos merőleges						GOOD
			{ { sin(incl), 0., cos(incl) }, { cos(incl), 0., -sin(incl) } },
			{ { cos(incl), 0., -sin(incl) }, { sin(incl), 0., cos(incl) } }, };
	double t, p, s;
	t = p = s = 0.;
	double fp = 0.5 * (1 + t * t) * cos(p) * cos(s) - t * sin(p) * sin(s);
	double fc = 0.5 * (1 + t * t) * cos(p) * sin(s) + t * sin(p) * cos(s);
short j;
detector_Struct max_Det;
	for (j = 0; j < 2; j++) {
		injParams.spin1x = chi[0] * where[j][0][0];
		injParams.spin1y = chi[0] * where[j][0][1];
		injParams.spin1z = chi[0] * where[j][0][2];
		injParams.spin2x = chi[1] * where[j][1][0];
		injParams.spin2y = chi[1] * where[j][1][1];
		injParams.spin2z = chi[1] * where[j][1][2];
		memset(&status, 0, sizeof(LALStatus));
		memset(&thewaveform1, 0, sizeof(CoherentGW));
		memset(&thewaveform2, 0, sizeof(CoherentGW));
		sprintf(PNString, "SpinQuadTaylor%s%s", "twoPointFivePN", "ALL");
		LALSnprintf(injParams.waveform, LIGOMETA_WAVEFORM_MAX * sizeof(CHAR),
				PNString);
		LALGenerateInspiral(&status, &thewaveform1, &injParams, &ppnParams);
		if (status.statusCode) {
			fprintf(stderr, "LALSQTPNWaveformTest: error generating waveform\n");
			return status.statusCode;
		}
		sprintf(PNString, "SpinQuadTaylor%s%s", "twoPointFivePN", "SS");
		LALSnprintf(injParams.waveform, LIGOMETA_WAVEFORM_MAX * sizeof(CHAR),
				PNString);
		//interface(&mystatus, &thewaveform, &injParams, &ppnParams);
		LALGenerateInspiral(&status, &thewaveform2, &injParams, &ppnParams);
		if (status.statusCode) {
			fprintf(stderr, "LALSQTPNWaveformTest: error generating waveform\n");
			return status.statusCode;
		}

		size_t max_Length = thewaveform1.f->data->length
				> thewaveform2.f->data->length ? thewaveform1.f->data->length
				: thewaveform2.f->data->length;
		freq_Max
				= thewaveform1.f->data->length > thewaveform2.f->data->length ? thewaveform1.f->data->data[thewaveform1.f->data->length
						- 1]
						: thewaveform2.f->data->data[thewaveform2.f->data->length
								- 1];
		static RandomInspiralSignalIn max_RandIn;
		max_RandIn.psd.length = max_Length;
		double max_Df = 1. / dt / max_RandIn.psd.length;
		max_RandIn.psd.data = (REAL8*) LALMalloc(sizeof(REAL8)
				* max_RandIn.psd.length);
		LALNoiseSpectralDensity(&status, &max_RandIn.psd, &LALLIGOIPsd, max_Df);
		double *max_PSD = fftw_malloc(max_Length * sizeof(double));
		for (i = 0; i < max_Length; i++) {
			max_PSD[i] = max_RandIn.psd.data[i];
		}
		double hp, hc, a1, a2, phi, shift;
		multi_Malloc(max_Length, &max_Det);
		memset(max_Det.t, 0, max_Det.length * sizeof(double));
		memset(max_Det.s, 0, max_Det.length * sizeof(double));
		for (i = 0; i < thewaveform1.f->data->length; i++) {
			a1 = thewaveform1.a->data->data[2 * i];
			a2 = thewaveform1.a->data->data[2 * i + 1];
			phi = thewaveform1.phi->data->data[i]
					- thewaveform1.phi->data->data[0];
			shift = thewaveform1.shift->data->data[i];
			hp = a1 * cos(shift) * cos(phi) - a2 * sin(shift) * sin(phi);
			hc = a1 * sin(shift) * cos(phi) + a2 * cos(shift) * sin(phi);
			max_Det.t[i] = fp * hp + fc * hc;
		}
		for (i = 0; i < thewaveform2.f->data->length; i++) {
			a1 = thewaveform1.a->data->data[2 * i];
			a2 = thewaveform1.a->data->data[2 * i + 1];
			phi = thewaveform1.phi->data->data[i]
					- thewaveform1.phi->data->data[0];
			shift = thewaveform1.shift->data->data[i];
			hp = a1 * cos(shift) * cos(phi) - a2 * sin(shift) * sin(phi);
			hc = a1 * sin(shift) * cos(phi) + a2 * cos(shift) * sin(phi);
			max_Det.s[i] = fp * hp + fc * hc;
//			printf("% -15.5lg % -15.5lg % -15.5lg % -15.5lg % -15.5lg\n", max_Det.s[i], fp, hp, fc, hc);
		}
		fftw_execute(max_Det.pt);
		fftw_execute(max_Det.ps);
		for (i = 0; i < 10/*thewaveform2.f->data->length*/; i++) {
			printf("X% -15.5lg % -15.5lg % -15.5lg % -15.5lg % -15.5lg\n", max_Det.s[i], max_Det.cs[0][i], max_Det.cs[1][i], max_Det.ct[0][i], max_Det.ct[1][i]);
		}
		XLALSQTPNDestroyCoherentGW(&thewaveform1);
		XLALSQTPNDestroyCoherentGW(&thewaveform2);
		double freq_Step, fr = 0.;
		freq_Step = 1. / (dt * max_Length);
		long max_Init_Fr = 0, max_Final_Fr = 0;
		while (fr < freq_Min) {
			fr += freq_Step;
			max_Final_Fr = ++max_Init_Fr;
		}
		while (fr < freq_Max) {
			fr += freq_Step;
			max_Final_Fr++;
		}
		double max_ts = scalar_freq(max_Det.ct, max_Det.cs, max_PSD,
				max_Init_Fr, max_Final_Fr);
		double max_tt = scalar_freq(max_Det.ct, max_Det.ct, max_PSD,
				max_Init_Fr, max_Final_Fr);
		double max_ss = scalar_freq(max_Det.cs, max_Det.cs, max_PSD,
				max_Init_Fr, max_Final_Fr);
		printf("filling: % -15.5lg % -15.5lg % -15.5lg: % -15.5lg\n", max_ts, max_tt, max_ss, max_ts / sqrt(max_tt * max_ss));
		fftw_free(max_PSD);
		multi_Free(&max_Det);
	}
	puts("Done.");
	return 0;
}
