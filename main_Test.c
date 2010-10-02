/** @file main_Test.c
 * @date 09/30/2010 08:38:34 PM
 * @author László Veréb
 * @brief
 */

#include "generator.h"

#include <lal/GenerateInspiral.h>
#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/LALSQTPNWaveformInterface.h>

NRCSID(LALSQTPNWAVEFORMTESTC, "$Id: LALSQTPNWaveformTest.c,v 0.1 2010/05/21");

int lalDebugLevel = 0; ///< the debug level

#define PREC "% -15.9lg "
#define PREC2 "% -15.9lg| "

#define DIR "test_Dir/"
int test(void);

int generate(long length);

int SQT_diff_ST(long length);

int main(int argc, char *argv[]) {
	SQT_diff_ST(atol(argv[1]));
	//	generate(atoi(argv[1]));
	puts("Done.");
	return 0;
}

void destroySTWave(CoherentGW waveform) {
	XLALDestroyREAL4Vector(waveform.f->data);
	LALFree(waveform.f);
	XLALDestroyREAL4Vector(waveform.shift->data);
	LALFree(waveform.shift);
	XLALDestroyREAL8Vector(waveform.phi->data);
	LALFree(waveform.phi);
	XLALDestroyREAL4VectorSequence(waveform.a->data);
	XLALFree(waveform.a);
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
	fclose(file_Gen);
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
		file_Gen = fopen(DIR"params.data", "a");
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
		fprintf(file_Gen, PREC2 PREC2 PREC2 PREC2 PREC2, act[i].dist, act[i].incl,
				act[i].F.dec, act[i].F.pol, act[i].F.phi);
		fflush(file_Gen);
		fprintf(file_Gen, "\n");
		fflush(file_Gen);
		fclose(file_Gen);
		injParams[0].mass1 = act[i].bh[0].m;
		injParams[0].mass2 = act[i].bh[1].m;
		injParams[0].spin1x = act[i].bh[0].chi[0];
		injParams[0].spin1y = act[i].bh[0].chi[1];
		injParams[0].spin1z = act[i].bh[0].chi[2];
		injParams[0].spin2x = act[i].bh[1].chi[0];
		injParams[0].spin2y = act[i].bh[1].chi[1];
		injParams[0].spin2z = act[i].bh[1].chi[2];
		injParams[0].qmParameter1 = 1.;
		injParams[0].qmParameter2 = 1.;
		injParams[0].inclination = act[i].incl;
		injParams[0].f_lower = 50.;
		injParams[0].f_final = 0.;
		injParams[0].distance = act[i].dist;
		injParams[0].polarization = 0;
		ppnParams.deltaT = 6.1035156e-05;
		memcpy(&injParams[1], &injParams[0], sizeof(SimInspiralTable));
		LALSnprintf(injParams[0].waveform,
				LIGOMETA_WAVEFORM_MAX * sizeof(CHAR),
				"SpinQuadTaylortwoPointFivePNSS");
		LALSnprintf(injParams[1].waveform,
				LIGOMETA_WAVEFORM_MAX * sizeof(CHAR),
				"SpinTaylortwoPointFivePN");
		for (j = 0; j < 2; j++) {
			memset(&status, 0, sizeof(LALStatus));
			memset(&waveform[j], 0, sizeof(CoherentGW));
			LALGenerateInspiral(&status, &waveform[j], &injParams[j],
					&ppnParams);
			if (status.statusCode) {
				fprintf(
						stderr,
						"%d: LALS(Q)TPNWaveformTest: error generating waveform\n",
						j);
			}
		}
		sprintf(file_Out_Name, DIR"gen%04ld.data", i);
		file_Out = fopen(file_Out_Name, "w");
		long shorter = waveform[0].f->data->length
				< waveform[1].f->data->length ? 0 : 1;
		for (k = 0; k < waveform[shorter].f->data->length; k++) {
			for (j = 0; j < 2; j++) {
				cshift = cos(waveform[j].shift->data->data[k]);
				sshift = sin(waveform[j].shift->data->data[k]);
				cphi = cos(waveform[j].phi->data->data[k]
						- waveform[j].phi->data->data[0]);
				sphi = sin(waveform[j].phi->data->data[k]
						- waveform[j].phi->data->data[0]);
				hp = waveform[j].a->data->data[2 * k] * cshift * cphi
						- waveform[j].a->data->data[2 * k + 1] * sshift * sphi;
				hc = waveform[j].a->data->data[2 * k] * sshift * cphi
						+ waveform[j].a->data->data[2 * k + 1] * cshift * sphi;
				act[i].F.F[0] = act[i].F.F[1] = sqrt(2.) / 2.;
				h[j] = act[i].F.F[0] * hp + act[i].F.F[1] * hc;
			}
			fprintf(file_Out, PREC PREC PREC PREC PREC"\n", k * ppnParams.deltaT, h[0], h[1], fabs(h[0]-h[1]), fabs(h[0]-h[1]) / h[0]);
		}
		XLALSQTPNDestroyCoherentGW(&waveform[0]);
		destroySTWave(waveform[1]);
		fclose(file_Out);
		//// t[i] = time() - t_I;
	}
	//	fclose(file_Gen);
	return 0.;
}

int generate(long length) {
	long i;
	binary_System sys, min, max;
	min.bh[0].m = min.bh[1].m = 3.;
	max.bh[0].m = 2 * (max.bh[1].m = 15.);
	min.M = 6.;
	max.M = 40.;
	min.eta = 0.01;
	max.eta = 0.25;
	min.bh[0].chi_Amp = min.bh[1].chi_Amp = 0.5;
	max.bh[0].chi_Amp = max.bh[1].chi_Amp = 1.;
	min.bh[0].cth = min.bh[1].cth = -(max.bh[0].cth = max.bh[1].cth = 1.);
	min.bh[0].phi = min.bh[1].phi = 0.;
	max.bh[0].phi = max.bh[1].phi = 2. * M_PI;
	max.dist = 10. * (min.dist = 1.);
	min.incl = 0.;
	max.incl = 2. * M_PI;
	min.F.dec = min.F.pol = min.F.phi = 0.;
	max.F.dec = max.F.pol = max.F.phi = 1.;
	srand(10);
	FILE *file = fopen(DIR"params.data", "w");
	for (i = 0; i < length; i++) {
		gen_Parameters(&sys, &min, &max, ETAM);
		fprintf(file, PREC PREC PREC PREC PREC, DBL_MIN, sys.M, sys.eta,
				sys.bh[0].m, sys.bh[1].m);
		fprintf(file, PREC PREC PREC, sys.bh[0].chi_Amp, sys.bh[0].cth,
				sys.bh[0].phi);
		fprintf(file, PREC PREC PREC, sys.bh[1].chi_Amp, sys.bh[1].cth,
				sys.bh[1].phi);
		fprintf(file, PREC PREC PREC PREC PREC, sys.dist, sys.incl, sys.F.dec,
				sys.F.pol, sys.F.phi);
		fprintf(file, "\n");
		fflush(file);
	}
	fclose(file);
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
	convert_etaM_m1m2(&sys, ETAM_TO_M1M2);
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
	char *filename[2] = { "2PNALL-47.out", "2PNSS-47.out" };
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
