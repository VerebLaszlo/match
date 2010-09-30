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

int main() {
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
