/** @file main_Test.c
 * @date 09/30/2010 08:38:34 PM
 * @author László Veréb
 * @brief
 */

#include "generator.h"
#include "match.h"
#include "detector.h"

#include <lal/LALNoiseModelsInspiral.h>
#include <lal/GenerateInspiral.h>
#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/LALSQTPNWaveformInterface.h>

NRCSID(LALSQTPNWAVEFORMTESTC, "$Id: LALSQTPNWaveformTest.c,v 0.1 2010/05/21");

int lalDebugLevel = 0; ///< the debug level

#define PREC2 "% -15.9lg| "

#define DIR "test_Dir/"
int test(void);

int generate(long length);

int calc_Match(long length);

int SQT_diff_ST(long length);

int measure_Time(long length);

int main(int argc, char *argv[]) {
	int x = argc;
	x = argc;
	calc_Match(atol(argv[1]));
	//measure_Time(atol(argv[1]));
	//SQT_diff_ST(atol(argv[1]));
	//generate(atoi(argv[1]));
	puts("Done.");
	return 0;
}

void destroySTWave(CoherentGW waveform) {
	if (waveform.f->data)
		XLALDestroyREAL4Vector(waveform.f->data);
	if (waveform.f)
		LALFree(waveform.f);
	if (waveform.shift->data)
		XLALDestroyREAL4Vector(waveform.shift->data);
	if (waveform.shift)
		LALFree(waveform.shift);
	if (waveform.phi->data)
		XLALDestroyREAL8Vector(waveform.phi->data);
	if (waveform.phi)
		LALFree(waveform.phi);
	if (waveform.a->data)
		XLALDestroyREAL4VectorSequence(waveform.a->data);
	if (waveform.a)
		XLALFree(waveform.a);
}

#define MODE "SpinTaylorthreePointFivePN"
#define PN "threePointFivePN"

int calc_Match(long length) {
	static LALStatus status;
	CoherentGW waveform[2];
	SimInspiralTable injParams[2];
	PPNParamStruc ppnParams;
	binary_System act;
	binary_System min, max;
	signalStruct signal, typ1, typ2, simp;
	static RandomInspiralSignalIn randIn;
	double best, worst, bestT, worstT;
	//// double t[length + 1];
	//// double t_I = time();
	/*FILE *file_Gen = fopen(DIR"params.data", "w");
	 fprintf(file_Gen, "#%-14s| %-15s| %-15s| %-15s| ", "M", "eta", "m_1", "m_2");
	 fprintf(file_Gen, "%-15s| %-15s| %-15s| %-15s| %-15s| %-15s| ", "chi_1",
	 "cth_1", "phi_1", "chi_1x", "chi_1y", "chi_1z");
	 fprintf(file_Gen, "%-15s| %-15s| %-15s| %-15s| %-15s| %-15s| ", "chi_2",
	 "cth_2", "phi_2", "chi_2x", "chi_2y", "chi_2z");
	 fprintf(file_Gen, "%-15s| %-15s| %-15s| %-15s| %-15s|\n", "dist", "incs",
	 "dec", "pol", "phi");*/
	//	fclose(file_Gen);
	//	FILE *file_Out;
	//	char file_Out_Name[50];
	memset(&injParams, 0, 2 * sizeof(SimInspiralTable));
	memset(&ppnParams, 0, sizeof(PPNParamStruc));
	min.bh[0].m = min.bh[1].m = min.bh[0].cth = min.bh[1].cth = min.bh[0].phi
			= min.bh[1].phi = min.bh[0].chi_Amp = min.bh[1].chi_Amp
					= min.bh[0].chi[0] = min.bh[0].chi[1] = min.bh[0].chi[2]
							= min.bh[1].chi[0] = min.bh[1].chi[1]
									= min.bh[1].chi[2] = min.M = min.eta
											= min.incl = min.dist = min.F.dec
													= min.F.pol = min.F.phi
															= -1000.;
	max.bh[0].m = max.bh[1].m = max.bh[0].cth = max.bh[1].cth = max.bh[0].phi
			= max.bh[1].phi = max.bh[0].chi_Amp = max.bh[1].chi_Amp
					= max.bh[0].chi[0] = max.bh[0].chi[1] = max.bh[0].chi[2]
							= max.bh[1].chi[0] = max.bh[1].chi[1]
									= max.bh[1].chi[2] = max.M = max.eta
											= max.incl = max.dist = max.F.dec
													= max.F.pol = max.F.phi
															= 1000.;
	check_Borders(&min, &max);

	double hp, hc, cphi, sphi, cshift, sshift;
	long i;
	unsigned long k;
	unsigned short j;
	//// t[0] = time() - t_I;
	for (i = 0; i < length; i++) {
		printf("%ld %ld\n", length, i);
		gen_Parameters(&act, &min, &max, ETAM);
		calc_Response_For_Detector(LH, act.F.dec, act.F.phi, act.F.pol,
				&act.F.F[0], &act.F.F[1]);
		//		file_Gen = fopen(DIR"params.data", "a");
		/*		fprintf(file_Gen, PREC2 PREC2 PREC2 PREC2, act.M, act.eta,
		 act.bh[0].m, act.bh[1].m);
		 fflush(file_Gen);
		 fprintf(file_Gen, PREC2 PREC2 PREC2, act.bh[0].chi_Amp,
		 act.bh[0].cth, act.bh[0].phi);
		 fflush(file_Gen);
		 fprintf(file_Gen, PREC2 PREC2 PREC2, act.bh[0].chi[0],
		 act.bh[0].chi[1], act.bh[0].chi[2]);
		 fflush(file_Gen);
		 fprintf(file_Gen, PREC2 PREC2 PREC2, act.bh[1].chi_Amp,
		 act.bh[1].cth, act.bh[1].phi);
		 fflush(file_Gen);
		 fprintf(file_Gen, PREC2 PREC2 PREC2, act.bh[1].chi[0],
		 act.bh[1].chi[1], act.bh[1].chi[2]);
		 fflush(file_Gen);
		 fprintf(file_Gen, PREC2 PREC2 PREC2 PREC2 PREC2, act.dist,
		 act.incl, act.F.dec, act.F.pol, act.F.phi);
		 fflush(file_Gen);
		 fprintf(file_Gen, "\n");
		 fflush(file_Gen);*/
		//		fclose(file_Gen);
		injParams[0].mass1 = act.bh[0].m = 34.;
		injParams[0].mass2 = act.bh[1].m = 3.4;
		injParams[0].spin1x = act.bh[0].chi[0];// = 0.0790765695;
		injParams[0].spin1y = act.bh[0].chi[1];// = 0.1589679868;
		injParams[0].spin1z = act.bh[0].chi[2];// = 0.9820794649;
		injParams[0].spin2x = act.bh[1].chi[0];// = 0.118388124;
		injParams[0].spin2y = act.bh[1].chi[1];// = -0.1715059157;
		injParams[0].spin2z = act.bh[1].chi[2];// = 0.9759989616;
		injParams[0].qmParameter1 = 1.;
		injParams[0].qmParameter2 = 1.;
		injParams[0].inclination = act.incl;// = 1.43;
		double freq_Min = injParams[0].f_lower = 40.;
		injParams[0].f_final = 0.;
		injParams[0].distance = act.dist;// = 1.;
		injParams[0].polarization = 0;
		ppnParams.deltaT = 6.1035156e-05;
		//		act.F.F[0] = act.F.F[1] = sqrt(2.) / 2.;
		//act.F.F[0] = act.F.F[1];
		memcpy(&injParams[1], &injParams[0], sizeof(SimInspiralTable));
		///***   PROBA   ***///
		//injParams[1].mass1 = 2. * injParams[0].mass1;	// valami problémát okoz
		///***   PROBA   ***///
		LALSnprintf(injParams[0].waveform,
				LIGOMETA_WAVEFORM_MAX * sizeof(CHAR), "SpinQuadTaylor"PN"SS");
		//LIGOMETA_WAVEFORM_MAX * sizeof(CHAR), "SpinQuadTaylor"PN"SELF");
		LALSnprintf(injParams[1].waveform,
				LIGOMETA_WAVEFORM_MAX * sizeof(CHAR), "SpinQuadTaylor"PN"ALL");
		//LIGOMETA_WAVEFORM_MAX * sizeof(CHAR), "SpinQuadTaylor"PN"SELF");
		memset(&status, 0, sizeof(LALStatus));
		memset(&waveform[0], 0, sizeof(CoherentGW));
		LALGenerateInspiral(&status, &waveform[0], &injParams[0], &ppnParams);
		if (status.statusCode) {
			fprintf(stderr, "LALSQTPNWaveformTest: error generating waveform\n");
			continue;
		}
		memset(&status, 0, sizeof(LALStatus));
		memset(&waveform[1], 0, sizeof(CoherentGW));
		LALGenerateInspiral(&status, &waveform[1], &injParams[1], &ppnParams);
		if (status.statusCode) {
			fprintf(stderr, "LALSTPNWaveformTest: error generating waveform\n");
			XLALSQTPNDestroyCoherentGW(&waveform[0]);
			continue;
		}
		//		sprintf(file_Out_Name, DIR"gen%04ld.data", i);
		//	file_Out = fopen(file_Out_Name, "w");

		long shorter = waveform[0].f->data->length
				< waveform[1].f->data->length ? 0 : 1;
		long length = waveform[!shorter].f->data->length;
		waveform[!shorter].f->data->length = length;
		create_Signal_Struct(&signal, waveform[!shorter].f->data->length);
		create_Signal_Struct(&typ1, waveform[!shorter].f->data->length);
		create_Signal_Struct(&typ2, waveform[!shorter].f->data->length);
		create_Signal_Struct(&simp, waveform[!shorter].f->data->length);
		randIn.psd.length = waveform[!shorter].f->data->length;
		double df = 1. / ppnParams.deltaT / randIn.psd.length;
		randIn.psd.data = (REAL8*) LALMalloc(sizeof(REAL8) * randIn.psd.length);
		LALNoiseSpectralDensity(&status, &randIn.psd, &LALLIGOIPsd, df);
		for (j = 0; j < randIn.psd.length; j++) {
			typ1.psd[j] = typ2.psd[j] = simp.psd[j] = signal.psd[j]
					= randIn.psd.data[j];
		}
		for (j = 0; j < 2; j++) {
			for (k = 0; k < waveform[j].f->data->length; k++) {
				cshift = cos(waveform[j].shift->data->data[k]);
				sshift = sin(waveform[j].shift->data->data[k]);
				cphi = cos(waveform[j].phi->data->data[k]);
				sphi = sin(waveform[j].phi->data->data[k]);
				hp = waveform[j].a->data->data[2 * k] * cshift * cphi
						- waveform[j].a->data->data[2 * k + 1] * sshift * sphi;
				hc = waveform[j].a->data->data[2 * k] * sshift * cphi
						+ waveform[j].a->data->data[2 * k + 1] * cshift * sphi;
				//				act.F.F[0] = act.F.F[1] = sqrt(2.) / 2.;
				//				signal.signal[2 * j][k] = act.F.F[0] * hp[0];
				//				signal.signal[2 * j + 1][k] = act.F.F[1] * hc[0];
				typ2.signal[2 * j][k] = typ1.signal[2 * j][k] = simp.signal[2
						* j][k] = signal.signal[2 * j][k] = hp;
				//act.F.F[0] * hp[0] + act.F.F[1] * hc[0];
				typ2.signal[2 * j + 1][k] = typ1.signal[2 * j + 1][k]
						= simp.signal[2 * j + 1][k]
								= signal.signal[2 * j + 1][k] = hc;
				//act.F.F[0] * hp[1] + act.F.F[1] * hc[1];
			}
			XLALSQTPNDestroyCoherentGW(&waveform[j]);
		}
		double freq_Max = (injParams[0].f_final + injParams[1].f_final) / 2.;
		double freq_Step, fr = 0.;
		freq_Step = 1. / (ppnParams.deltaT * randIn.psd.length);
		long minfr = 0, maxfr = 0;
		while (fr < freq_Min) {
			fr += freq_Step;
			maxfr = ++minfr;
		}
		while (fr < freq_Max) {
			fr += freq_Step;
			maxfr++;
		}
		//		proba(&signal, minfr, maxfr);
		//		printf("Simp = "PREC"\n", match_Simple(&simp, minfr, maxfr));fflush(stdout);
		printf("T  = "PREC"\n", match_Typical(&typ1, minfr, maxfr, NONE));
		fflush(stdout);
		printf("TT = "PREC"\n", match_Typical(&typ2, minfr, maxfr, TIME));
		fflush(stdout);
		calc_Overlap(&best, &worst, &signal, minfr, maxfr);
		printf("B = "PREC"\nW = "PREC"\n", best, worst);
		calc_Overlap_Time(&bestT, &worstT, &simp, minfr, maxfr);
		printf("BT = "PREC"\nWT = "PREC"\n", bestT, worstT);
		destroy_Signal_Struct(&simp);
		destroy_Signal_Struct(&typ1);
		destroy_Signal_Struct(&typ2);
		destroy_Signal_Struct(&signal);
		//destroySTWave(waveform[1]);
		XLALFree(randIn.psd.data);
		//		fclose(file_Out);
		//// t[i] = time() - t_I;
	}
	//	fclose(file_Gen);
	return 0.;
}

#include <time.h>

int measure_Time(long length) {
	double time[length / 100 + 10];
	static LALStatus status;
	CoherentGW waveform;
	SimInspiralTable injParams;
	PPNParamStruc ppnParams;
	binary_System act, min, max;
	memset(&injParams, 0, sizeof(SimInspiralTable));
	memset(&ppnParams, 0, sizeof(PPNParamStruc));
	min.bh[0].m = min.bh[1].m = min.bh[0].cth = min.bh[1].cth = min.bh[0].phi
			= min.bh[1].phi = min.bh[0].chi_Amp = min.bh[1].chi_Amp
					= min.bh[0].chi[0] = min.bh[0].chi[1] = min.bh[0].chi[2]
							= min.bh[1].chi[0] = min.bh[1].chi[1]
									= min.bh[1].chi[2] = min.M = min.eta
											= min.incl = min.dist = min.F.dec
													= min.F.pol = min.F.phi
															= -1000.;
	max.bh[0].m = max.bh[1].m = max.bh[0].cth = max.bh[1].cth = max.bh[0].phi
			= max.bh[1].phi = max.bh[0].chi_Amp = max.bh[1].chi_Amp
					= max.bh[0].chi[0] = max.bh[0].chi[1] = max.bh[0].chi[2]
							= max.bh[1].chi[0] = max.bh[1].chi[1]
									= max.bh[1].chi[2] = max.M = max.eta
											= max.incl = max.dist = max.F.dec
													= max.F.pol = max.F.phi
															= 1000.;
	check_Borders(&min, &max);
	min.M = 15.;
	long i;
	clock_t t_I = clock();
	for (i = 0; i < length; i++) {
		gen_Parameters(&act, &min, &max, ETAM);
		injParams.mass1 = act.bh[0].m;
		injParams.mass2 = act.bh[1].m;
		injParams.spin1x = act.bh[0].chi[0];
		injParams.spin1y = act.bh[0].chi[1];
		injParams.spin1z = act.bh[0].chi[2];
		injParams.spin2x = act.bh[1].chi[0];
		injParams.spin2y = act.bh[1].chi[1];
		injParams.spin2z = act.bh[1].chi[2];
		injParams.qmParameter1 = 1.;
		injParams.qmParameter2 = 1.;
		injParams.inclination = act.incl;
		injParams.f_lower = 40.;
		injParams.f_final = 0.;
		injParams.distance = act.dist;
		injParams.polarization = 0;
		ppnParams.deltaT = 6.1035156e-05;
		LALSnprintf(injParams.waveform, LIGOMETA_WAVEFORM_MAX * sizeof(CHAR),
				MODE);
		memset(&status, 0, sizeof(LALStatus));
		memset(&waveform, 0, sizeof(CoherentGW));
		LALGenerateInspiral(&status, &waveform, &injParams, &ppnParams);
		if (status.statusCode) {
			fprintf(stderr, MODE": error generating waveform\n");
			continue;
		}
		XLALSQTPNDestroyCoherentGW(&waveform);
		//destroySTWave(waveform);
		if (i % 100 == 99) {
			time[i / 100] = (double) (clock() - t_I) / CLOCKS_PER_SEC;
			printf("%ld: %lg\n", i + 1, time[i / 100]);
		}
	}
	FILE *file = fopen(DIR"time35PNST.data", "w");
	for (i = 0; i < length / 100; i++) {
		fprintf(file, PREC"\n", time[i]);
		fflush(stdout);
	}
	fclose(file);
	return 0;
}

int SQT_diff_ST(long length) {
	static LALStatus status;
	CoherentGW waveform[2];
	SimInspiralTable injParams[2];
	PPNParamStruc ppnParams;
	binary_System act[length];
	binary_System min, max;
	//// double t[length + 1];
	//// double t_I = time();
	FILE *file_Gen = fopen(DIR"params.data", "w");
	fprintf(file_Gen, "#%-14s| %-15s| %-15s| %-15s| ", "M", "eta", "m_1", "m_2");
	fprintf(file_Gen, "%-15s| %-15s| %-15s| %-15s| %-15s| %-15s| ", "chi_1",
			"cth_1", "phi_1", "chi_1x", "chi_1y", "chi_1z");
	fprintf(file_Gen, "%-15s| %-15s| %-15s| %-15s| %-15s| %-15s| ", "chi_2",
			"cth_2", "phi_2", "chi_2x", "chi_2y", "chi_2z");
	fprintf(file_Gen, "%-15s| %-15s| %-15s| %-15s| %-15s|\n", "dist", "incs",
			"dec", "pol", "phi");
	//	fclose(file_Gen);
	FILE *file_Out;
	char file_Out_Name[50];
	memset(&injParams, 0, 2 * sizeof(SimInspiralTable));
	memset(&ppnParams, 0, sizeof(PPNParamStruc));
	min.bh[0].m = min.bh[1].m = min.bh[0].cth = min.bh[1].cth = min.bh[0].phi
			= min.bh[1].phi = min.bh[0].chi_Amp = min.bh[1].chi_Amp
					= min.bh[0].chi[0] = min.bh[0].chi[1] = min.bh[0].chi[2]
							= min.bh[1].chi[0] = min.bh[1].chi[1]
									= min.bh[1].chi[2] = min.M = min.eta
											= min.incl = min.dist = min.F.dec
													= min.F.pol = min.F.phi
															= -1000.;
	max.bh[0].m = max.bh[1].m = max.bh[0].cth = max.bh[1].cth = max.bh[0].phi
			= max.bh[1].phi = max.bh[0].chi_Amp = max.bh[1].chi_Amp
					= max.bh[0].chi[0] = max.bh[0].chi[1] = max.bh[0].chi[2]
							= max.bh[1].chi[0] = max.bh[1].chi[1]
									= max.bh[1].chi[2] = max.M = max.eta
											= max.incl = max.dist = max.F.dec
													= max.F.pol = max.F.phi
															= 1000.;
	check_Borders(&min, &max);
	double h[2], hp, hc, cphi, sphi, cshift, sshift;
	long i;
	unsigned long k;
	short j;
	//// t[0] = time() - t_I;
	for (i = 0; i < length; i++) {
		printf("%ld %ld\n", length, i);
		gen_Parameters(&act[i], &min, &max, ETAM);
		calc_Response_For_Detector(LH, act[i].F.dec, act[i].F.phi,
				act[i].F.pol, &act[i].F.F[0], &act[i].F.F[1]);
		//		file_Gen = fopen(DIR"params.data", "a");
		fprintf(file_Gen, PREC2 PREC2 PREC2 PREC2, act[i].M, act[i].eta,
				act[i].bh[0].m, act[i].bh[1].m);
		fflush(file_Gen);
		fprintf(file_Gen, PREC2 PREC2 PREC2, act[i].bh[0].chi_Amp,
				act[i].bh[0].cth, act[i].bh[0].phi);
		fflush(file_Gen);
		fprintf(file_Gen, PREC2 PREC2 PREC2, act[i].bh[0].chi[0],
				act[i].bh[0].chi[1], act[i].bh[0].chi[2]);
		fflush(file_Gen);
		fprintf(file_Gen, PREC2 PREC2 PREC2, act[i].bh[1].chi_Amp,
				act[i].bh[1].cth, act[i].bh[1].phi);
		fflush(file_Gen);
		fprintf(file_Gen, PREC2 PREC2 PREC2, act[i].bh[1].chi[0],
				act[i].bh[1].chi[1], act[i].bh[1].chi[2]);
		fflush(file_Gen);
		fprintf(file_Gen, PREC2 PREC2 PREC2 PREC2 PREC2, act[i].dist,
				act[i].incl, act[i].F.dec, act[i].F.pol, act[i].F.phi);
		fflush(file_Gen);
		fprintf(file_Gen, "\n");
		fflush(file_Gen);
		//		fclose(file_Gen);
		act[i].eta = 0.24;
		act[i].M = 4.6 / pow(act[i].eta, 3. / 5.);
		injParams[0].mass1 = act[i].bh[0].m = 4.5;
		injParams[0].mass2 = act[i].bh[1].m = 4.5;
		act[i].bh[0].chi_Amp = 0.6;
		act[i].bh[0].cth = cos(0.);
		act[i].bh[0].phi = 80. * M_PI / 180.;
		act[i].bh[0].chi_Amp = 0.7;
		act[i].bh[0].cth = cos(160. * M_PI / 180.);
		act[i].bh[0].phi = 80. * M_PI / 180.;
		convert_Angles_Components(&act[i], ANGLE_TO_COMP);
		injParams[0].spin1x = act[i].bh[0].chi[0];
		injParams[0].spin1y = act[i].bh[0].chi[1];
		injParams[0].spin1z = act[i].bh[0].chi[2];
		injParams[0].spin2x = act[i].bh[1].chi[0];
		injParams[0].spin2y = act[i].bh[1].chi[1];
		injParams[0].spin2z = act[i].bh[1].chi[2];
		injParams[0].qmParameter1 = 1.;
		injParams[0].qmParameter2 = 1.;
		injParams[0].inclination = act[i].incl = 15. * M_PI / 180.;
		injParams[0].f_lower = 40.;
		injParams[0].f_final = 0.;
		injParams[0].distance = act[i].dist = 47.;
		injParams[0].polarization = 0;
		injParams[0].coa_phase = 150. * M_PI / 180.;
		ppnParams.deltaT = 6.1035156e-05;
		memcpy(&injParams[1], &injParams[0], sizeof(SimInspiralTable));
		LALSnprintf(injParams[0].waveform,
				LIGOMETA_WAVEFORM_MAX * sizeof(CHAR), "SpinQuadTaylor"PN"SS");
		LALSnprintf(injParams[1].waveform,
				LIGOMETA_WAVEFORM_MAX * sizeof(CHAR), "SpinQuadTaylor"PN"SS");
		memset(&status, 0, sizeof(LALStatus));
		memset(&waveform[0], 0, sizeof(CoherentGW));
		LALGenerateInspiral(&status, &waveform[0], &injParams[0], &ppnParams);
		if (status.statusCode) {
			fprintf(stderr, "LALSQTPNWaveformTest: error generating waveform\n");
			continue;
		}
		memset(&status, 0, sizeof(LALStatus));
		memset(&waveform[1], 0, sizeof(CoherentGW));
		LALGenerateInspiral(&status, &waveform[1], &injParams[1], &ppnParams);
		if (status.statusCode) {
			fprintf(stderr, "LALSTPNWaveformTest: error generating waveform\n");
			XLALSQTPNDestroyCoherentGW(&waveform[0]);
			continue;
		}
		sprintf(file_Out_Name, DIR"gen%04ld.data", i);
		file_Out = fopen(file_Out_Name, "w");
		long shorter = waveform[0].f->data->length
				< waveform[1].f->data->length ? 0 : 1;
		for (k = 0; k < waveform[shorter].f->data->length; k++) {
			for (j = 0; j < 1; j++) {
				cshift = cos(waveform[j].shift->data->data[k]);
				sshift = sin(waveform[j].shift->data->data[k]);
				cphi = cos(waveform[j].phi->data->data[k]);
				//						- waveform[j].phi->data->data[0]);
				sphi = sin(waveform[j].phi->data->data[k]);
				//						- waveform[j].phi->data->data[0]);
				hp = waveform[j].a->data->data[2 * k] * cshift * cphi
						- waveform[j].a->data->data[2 * k + 1] * sshift * sphi;
				hc = waveform[j].a->data->data[2 * k] * sshift * cphi
						+ waveform[j].a->data->data[2 * k + 1] * cshift * sphi;
				//				act[i].F.F[0] = act[i].F.F[1] = sqrt(2.) / 2.;
				//				h[j] = act[i].F.F[0] * hp + act[i].F.F[1] * hc;
			}
			fprintf(file_Out, PREC PREC PREC PREC PREC"\n", k
					* ppnParams.deltaT, hp, hc, hp - hc, (hp - hc) / hp);
		}
		for (; k < waveform[!shorter].f->data->length; k++) {
			cshift = cos(waveform[!shorter].shift->data->data[k]);
			sshift = sin(waveform[!shorter].shift->data->data[k]);
			cphi = cos(waveform[!shorter].phi->data->data[k]
					- waveform[!shorter].phi->data->data[0]);
			sphi = sin(waveform[!shorter].phi->data->data[k]
					- waveform[!shorter].phi->data->data[0]);
			hp = waveform[!shorter].a->data->data[2 * k] * cshift * cphi
					- waveform[!shorter].a->data->data[2 * k + 1] * sshift
							* sphi;
			hc = waveform[!shorter].a->data->data[2 * k] * sshift * cphi
					+ waveform[!shorter].a->data->data[2 * k + 1] * cshift
							* sphi;
			act[i].F.F[0] = act[i].F.F[1] = sqrt(2.) / 2.;
			h[0] = act[i].F.F[0] * hp + act[i].F.F[1] * hc;
			fprintf(file_Out, PREC PREC PREC PREC PREC"\n", k
					* ppnParams.deltaT, h[0], 0., fabs(h[0]), fabs(h[0]) / h[0]);
		}
		XLALSQTPNDestroyCoherentGW(&waveform[0]);
		destroySTWave(waveform[1]);
		fclose(file_Out);
		//// t[i] = time() - t_I;
	}
	fclose(file_Gen);
	return 0.;
}

int generate(long length) {
	static LALStatus status;
	CoherentGW waveform[2];
	SimInspiralTable injParams[2];
	PPNParamStruc ppnParams;
	signalStruct signal;
	binary_System act, min, max;
	memset(&injParams, 0, 2 * sizeof(SimInspiralTable));
	memset(&ppnParams, 0, sizeof(PPNParamStruc));
	///***   PROBA   ***///
	/* 	act.bh[0].chi[0] = 0.0790765695;
	 act.bh[0].chi[1] = 0.1589679868;
	 act.bh[0].chi[2] = 0.9820794649;
	 act.bh[1].chi[0] = 0.118388124;
	 act.bh[1].chi[1] = -0.1715059157;
	 act.bh[1].chi[2] = 0.9759989616;
	 convert_Angles_Components(&act, COMP_TO_ANGLE);
	 printf(PREC PREC PREC"\n", act.bh[0].chi_Amp, act.bh[0].phi, act.bh[0].cth);
	 printf(PREC PREC PREC"\n", act.bh[1].chi_Amp, act.bh[1].phi, act.bh[1].cth);
	 exit(-1);
	 *////***   PROBA   ***///
	min.bh[0].m = min.bh[1].m = min.bh[0].cth = min.bh[1].cth = min.bh[0].phi
			= min.bh[1].phi = min.bh[0].chi_Amp = min.bh[1].chi_Amp
					= min.bh[0].chi[0] = min.bh[0].chi[1] = min.bh[0].chi[2]
							= min.bh[1].chi[0] = min.bh[1].chi[1]
									= min.bh[1].chi[2] = min.M = min.eta
											= min.incl = min.dist = min.F.dec
													= min.F.pol = min.F.phi
															= -1000.;
	max.bh[0].m = max.bh[1].m = max.bh[0].cth = max.bh[1].cth = max.bh[0].phi
			= max.bh[1].phi = max.bh[0].chi_Amp = max.bh[1].chi_Amp
					= max.bh[0].chi[0] = max.bh[0].chi[1] = max.bh[0].chi[2]
							= max.bh[1].chi[0] = max.bh[1].chi[1]
									= max.bh[1].chi[2] = max.M = max.eta
											= max.incl = max.dist = max.F.dec
													= max.F.pol = max.F.phi
															= 1000.;
	check_Borders(&min, &max);
	FILE *file = NULL;
	FILE *gen = fopen("generating/params.data", "w");
	char file_Name[50];
	double hp, hc, cphi, sphi, cshift, sshift;
	long i;
	unsigned long k;
	short j, longest;
	srand(10);
	for (i = 0; i < length; i++) {
		printf("%ld / %ld\n", i + 1, length);
		fflush(stdout);
		sprintf(file_Name, "generating/%04ld.out", i);
		file = fopen(file_Name, "w");
		gen_Parameters(&act, &min, &max, ETAM);
		/*		act.dist = 46.5;
		 act.incl = 15. / 180. * M_PI;
		 act.F.dec = -46. / 180. * M_PI;
		 act.F.pol = 74. / 180. * M_PI;
		 act.F.phi = 7.9 * 60. * 60.;	// h -> sec
		 calc_Response_For_Detector(LH, act.F.dec, act.F.phi, act.F.pol,
		 &act.F.F[0], &act.F.F[1]);
		 act.bh[0].chi_Amp = 0.6;
		 act.bh[1].chi_Amp = 0.7;
		 //act.bh[0].phi = acos((.1) / (sin(act.incl) *));
		 act.M = 4.6 / pow(0.24, 3. / 5.);*/
		injParams[0].mass1 = act.bh[0].m = 34;//(1. + sqrt(1. - 4. * act.eta)) * act.M / 2.;
		injParams[0].mass2 = act.bh[1].m = 3.4;//(1. - sqrt(1. - 4. * act.eta)) * act.M / 2.;
		injParams[0].spin1x = act.bh[0].chi[0];
		injParams[0].spin1y = act.bh[0].chi[1];
		injParams[0].spin1z = act.bh[0].chi[2];
		injParams[0].spin2x = act.bh[1].chi[0];
		injParams[0].spin2y = act.bh[1].chi[1];
		injParams[0].spin2z = act.bh[1].chi[2];
		injParams[0].qmParameter1 = 1.;
		injParams[0].qmParameter2 = 1.;
		injParams[0].inclination = act.incl;
		injParams[0].f_lower = 40.;
		injParams[0].f_final = 0.;
		injParams[0].distance = act.dist;
		injParams[0].polarization = 0;
		print_Binary_System(&act, gen);
		ppnParams.deltaT = 6.1035156e-05;
		memcpy(&injParams[1], &injParams[0], sizeof(SimInspiralTable));
		injParams[0].coa_phase = M_PI / 2.;
		injParams[1].coa_phase = 0.;
		LALSnprintf(injParams[0].waveform,
				LIGOMETA_WAVEFORM_MAX * sizeof(CHAR), "SpinQuadTaylor"PN"NO");
		LALSnprintf(injParams[1].waveform,
				LIGOMETA_WAVEFORM_MAX * sizeof(CHAR), "SpinQuadTaylor"PN"NO");
		memset(&status, 0, sizeof(LALStatus));
		memset(&waveform[0], 0, sizeof(CoherentGW));
		LALGenerateInspiral(&status, &waveform[0], &injParams[0], &ppnParams);
		if (status.statusCode) {
			fprintf(stderr, "LALSQTPNWaveformTest: error generating waveform\n");
			continue;
		}
		memset(&status, 0, sizeof(LALStatus));
		memset(&waveform[1], 0, sizeof(CoherentGW));
		LALGenerateInspiral(&status, &waveform[1], &injParams[1], &ppnParams);
		if (status.statusCode) {
			fprintf(stderr, "LALSTPNWaveformTest: error generating waveform\n");
			XLALSQTPNDestroyCoherentGW(&waveform[0]);
			continue;
		}
		longest = waveform[0].f->data->length < waveform[1].f->data->length ? 0
				: 1;
		create_Signal_Struct(&signal, waveform[longest].f->data->length);
		for (j = 0; j < 2; j++) {
			for (k = 0; k < waveform[j].f->data->length; k++) {
				cshift = cos(waveform[j].shift->data->data[k]);
				sshift = sin(waveform[j].shift->data->data[k]);
				cphi = cos(waveform[j].phi->data->data[k]
						- waveform[j].phi->data->data[0]);// + j * M_PI / 2);
				sphi = sin(waveform[j].phi->data->data[k]
						- waveform[j].phi->data->data[0]);// + j * M_PI / 2);
				hp = waveform[j].a->data->data[2 * k] * cshift * cphi
						- waveform[j].a->data->data[2 * k + 1] * sshift * sphi;
				hc = waveform[j].a->data->data[2 * k] * sshift * cphi
						+ waveform[j].a->data->data[2 * k + 1] * cshift * sphi;
				signal.signal[2 * j][k] = act.F.F[0] * hp + act.F.F[1] * hc;
			}
			XLALSQTPNDestroyCoherentGW(&waveform[j]);
		}
		puts("R");
		for (k = 0; k < (unsigned long)signal.size; k++) {
			fprintf(file, PREC PREC PREC PREC"\n", k * ppnParams.deltaT,
					signal.signal[0][k], signal.signal[2][k],
					signal.signal[0][k] - signal.signal[2][k]);
			fflush(file);
		}
		destroy_Signal_Struct(&signal);
		fclose(file);
	}
	fclose(gen);
	return 0;
}

int test(void) {
	static LALStatus status;
	CoherentGW waveform[2];
	SimInspiralTable injParams[2];
	PPNParamStruc ppnParams;
	binary_System sys;
	memset(injParams, 0, 2 * sizeof(SimInspiralTable));
	memset(&ppnParams, 0, sizeof(PPNParamStruc));
	memset(&sys, 0, sizeof(binary_System));
	sys.M = 4.58;
	sys.eta = 0.242;
	sys.dist = 46.;
	sys.incl = DEG_TO_RAD(-47.);
	sys.bh[0].chi_Amp = 0.55;
	sys.bh[0].phi = DEG_TO_RAD(85.);
	sys.bh[0].cth = cos(DEG_TO_RAD(-60.));
	sys.bh[1].chi_Amp = 0.65;
	sys.bh[1].phi = DEG_TO_RAD(90.);
	sys.bh[1].cth = cos(DEG_TO_RAD(160.));
	convert_Angles_Components(&sys, ANGLE_TO_COMP);
	convert_etaM_m1m2(&sys, FROM_ETAM);
	injParams[0].mass1 = sys.bh[0].m;
	injParams[0].mass2 = sys.bh[1].m;
	injParams[0].spin1x = sys.bh[0].chi[0];
	injParams[0].spin1y = sys.bh[0].chi[1];
	injParams[0].spin1z = sys.bh[0].chi[2];
	injParams[0].spin2x = sys.bh[1].chi[0];
	injParams[0].spin2y = sys.bh[1].chi[1];
	injParams[0].spin2z = sys.bh[1].chi[2];
	injParams[0].qmParameter1 = 1.;
	injParams[0].qmParameter2 = 1.;
	injParams[0].inclination = sys.incl;
	injParams[0].f_lower = 50.;
	injParams[0].f_final = 0.;
	injParams[0].distance = sys.dist;
	injParams[0].polarization = 0;
	ppnParams.deltaT = 6.1035156e-05;
	memcpy(&injParams[1], &injParams[0], sizeof(SimInspiralTable));
	LALSnprintf(injParams[0].waveform, LIGOMETA_WAVEFORM_MAX * sizeof(CHAR),
			"SpinQuadTaylortwoPointFivePNALL");
	LALSnprintf(injParams[1].waveform, LIGOMETA_WAVEFORM_MAX * sizeof(CHAR),
			"SpinTaylortwoPointFivePN");
	char *filename[2] =
		{ "2PNALL-47.out", "2PNSS-47.out" };
	FILE *file;
	memset(&status, 0, sizeof(LALStatus));
	memset(waveform, 0, 2 * sizeof(CoherentGW));
	long j, length;
	double dt, a1, a2, phi, shift;
	short i;
	for (i = 0; i < 2; i++) {
		LALGenerateInspiral(&status, &waveform[i], &injParams[i], &ppnParams);
		if (status.statusCode) {
			fprintf(stderr,
					"%d: LALS(Q)TPNWaveformTest: error generating waveform\n",
					i);
			return status.statusCode;
		}
		file = fopen(filename[i], "w");
		length = waveform[i].f->data->length;
		dt = waveform[i].phi->deltaT;
		for (j = 0; j < length; j++) {
			a1 = waveform[i].a->data->data[2 * j];
			a2 = waveform[i].a->data->data[2 * j + 1];
			phi = waveform[i].phi->data->data[j]
					- waveform[i].phi->data->data[0];
			shift = waveform[i].shift->data->data[j];

			fprintf(file, PREC PREC PREC"\n", j * dt, a1 * cos(shift)
					* cos(phi) - a2 * sin(shift) * sin(phi), a1 * sin(shift)
					* cos(phi) + a2 * cos(shift) * sin(phi));
		}
		fclose(file);
		XLALSQTPNDestroyCoherentGW(&waveform[i]);
	}
	puts("Done.");
	LALCheckMemoryLeaks();
	return 0;
}
