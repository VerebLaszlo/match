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

#define PREC "% -15.10lg "

#define DIR "test_Dir/"
int test(void);

int generate(long count);

int main(int argc, char *argv[]) {
	generate(atoi(argv[1]));
	puts("Done.");
	return 0;
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
	min.incl = 0.; max.incl = 2. * M_PI;
	min.F.dec = min.F.pol = min.F.phi = 0.;
	max.F.dec = max.F.pol = max.F.phi = 1.;
	srand(10);
	FILE *file = fopen(DIR"params.data", "w");
	for (i = 0; i < length; i++) {
		gen_Parameters(&sys, &min, &max, ETAM);
		fprintf(file, PREC PREC PREC PREC PREC, DBL_MIN, sys.M, sys.eta, sys.bh[0].m, sys.bh[1].m);
		fprintf(file, PREC PREC PREC, sys.bh[0].chi_Amp, sys.bh[0].cth, sys.bh[0].phi);
		fprintf(file, PREC PREC PREC, sys.bh[1].chi_Amp, sys.bh[1].cth, sys.bh[1].phi);
		fprintf(file, PREC PREC PREC PREC PREC, sys.dist, sys.incl, sys.F.dec, sys.F.pol, sys.F.phi);
		fprintf(file, "\n");fflush(file);
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
	char *filename[2] = {"2PNALL-47.out", "2PNSS-47.out"};
	FILE *file;
	memset(&status, 0, sizeof(LALStatus));
	memset(waveform, 0, 2 * sizeof(CoherentGW));
	long j, length;
	double dt, a1, a2, phi, shift;
	short i;
	for (i = 0; i < 2; i++) {
		LALGenerateInspiral(&status, &waveform[i], &injParams[i], &ppnParams);
		if (status.statusCode) {
			fprintf(stderr, "%d: LALS(Q)TPNWaveformTest: error generating waveform\n", i);
			return status.statusCode;
		}
		file = fopen(filename[i], "w");
		length = waveform[i].f->data->length;
		dt      = waveform[i].phi->deltaT;
	    for(j = 0; j < length; j++) {
	        a1  = waveform[i].a->data->data[2*j];
	        a2  = waveform[i].a->data->data[2*j+1];
	        phi     = waveform[i].phi->data->data[j] - waveform[i].phi->data->data[0];
	        shift   = waveform[i].shift->data->data[j];

	        fprintf(file, PREC PREC PREC"\n",
	            j*dt,
	            a1*cos(shift)*cos(phi) - a2*sin(shift)*sin(phi),
	            a1*sin(shift)*cos(phi) + a2*cos(shift)*sin(phi));
	    }
		fclose(file);
		XLALSQTPNDestroyCoherentGW(&waveform[i]);
	}
	puts("Done.");
	LALCheckMemoryLeaks();
	return 0;
}
