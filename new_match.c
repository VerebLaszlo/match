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

double pi = LAL_PI;

int main(int argc, char *argv[]) {
	static LALStatus status;
	CoherentGW waveform[2];
	SimInspiralTable injParams;
	PPNParamStruc ppnParams;
	detector_Struct det;
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


	double freq_i, freq_f, dt;
	freq_i = injParams.f_lower;
	dt = ppnParams.deltaT;
	// antenna függvények kiszámolása
	double t, p, s, fp, fc;
	t = p = s = 0.;
	fp = 0.5 * (1 + t * t) * cos(p) * cos(s) - t * sin(p) * sin(s);
	fc = 0.5 * (1 + t * t) * cos(p) * sin(s) + t * sin(p) * cos(s);
	double hp, hc, a1, a2, phi, shift;
	double *data[2];
	long i, length;
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
		/*multi_Malloc(length, &det);
		data[0] = det.t;
		data[1] = det.s;*/
		for (second = 0; second < 2; second++) {
			memset(data[second], 0, det.length * sizeof(double));
			/*for (i = 0; i < waveform[second].f->data->length; i++) {
				a1 = waveform[second].a->data->data[2 * i];
				a2 = waveform[second].a->data->data[2 * i + 1];
				phi = waveform[second].phi->data->data[i] - waveform[second].phi->data->data[0];
				shift = waveform[second].shift->data->data[i];
				hp = a1 * cos(shift) * cos(phi) - a2 * sin(shift) * sin(phi);
				hc = a1 * sin(shift) * cos(phi) + a2 * cos(shift) * sin(phi);
				data[second][i] = fp * hp + fc * hc;
			}*/
			XLALSQTPNDestroyCoherentGW(&waveform[second]);
		}
		//multi_Free(&det);
	}
	puts("Done.");
	return 0;
}
