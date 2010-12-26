/*
 * main_Generator.c
 *
 *  Created on: Sep 29, 2010
 *      Author: vereb
 */

//#include "confuse-parser.h"

#define restrict __restrict__
#include "generator.h"
#include "match.h"
#include "detector.h"
#include "match_Multi.h"

#include <lal/LALNoiseModelsInspiral.h>
#include <lal/GenerateInspiral.h>
#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/LALSQTPNWaveformInterface.h>

double pi = M_PI;

double theta[2] =
	{ 3.075991792708861360949868, 2.316140092774243264273082 };

double varphi[2] =
	{ 1.371348580001392480909317, 3.849329592965536228632573 };

void mod_Spins(void);

void generate(void);

extern double maxT;
extern short maxTindex;

int main(int argc, char *argv[]) {
	generate();
	//mod_Spins();
	puts("Done!");
	return 0;
}

void mod_Spins(void) {
	char dir[50];
	sprintf(dir, "data_TDK/spinConf/");
	program_Params params;
	memset(&params, 0, sizeof(program_Params));
	params.freq_Sampling = 16380.;
	params.time_Sampling = 1. / params.freq_Sampling;
	params.freq_Initial = 40.;
	params.freq_Final = 0.;
	short i, j, k, l, index, step = 20;
	short big_Step = 1. * step;
	short length = big_Step * big_Step;
	double amp_Start = 0.998;
	double amp_Step = amp_Start / (double) step * 2.;
	double amp[2] =
		{ amp_Start - amp_Step, amp_Start - amp_Step };
	binary_System *act = malloc(length * sizeof(binary_System));
	memset(act, 0, length * sizeof(binary_System));
	for (i = 0; i < big_Step; i++) {
		amp[0] += amp_Step;
		for (j = 0; j < big_Step; j++) {
			amp[1] += amp_Step;
			//printf("% -20.15lg % -20.15lg\n", amp[0], amp[1]);
			index = i * big_Step + j;
			act[index].F.dec = 0.;
			act[index].F.pol = 0.;
			act[index].F.alpha = 0.;
			act[index].F.gmst = 0.;
			calc_Response_For_Detector(LH, &act[index]);
			act[index].incl = 1.34;
			act[index].dist = 1.;
			act[index].bh[0].m = 34.;
			act[index].bh[1].m = 3.4;
			convert_Masses(&act[index], FROM_M1M2);
			for (k = 0; k < 2; k++) {
				act[index].bh[k].chi_Amp = amp[k];
				act[index].bh[k].theta = theta[k];
				act[index].bh[k].varphi = varphi[k];
			}
			convert_Spins(&act[index], FROM_THETA_VPHI);
			/*printf("% -20.15lg % -20.15lg % -20.15lg % -20.15lg % -20.15lg\n",
			 act[index].bh[0].chi_Amp, sqrt(SQR(act[index].bh[0].chi[0])
			 + SQR(act[index].bh[0].chi[1])
			 + SQR(act[index].bh[0].chi[2])),
			 act[index].bh[0].chi[0], act[index].bh[0].chi[1],
			 act[index].bh[0].chi[2]);
			 printf("% -20.15lg % -20.15lg % -20.15lg % -20.15lg % -20.15lg\n",
			 act[index].bh[1].chi_Amp, sqrt(SQR(act[index].bh[1].chi[0])
			 + SQR(act[index].bh[1].chi[1])
			 + SQR(act[index].bh[1].chi[2])),
			 act[index].bh[1].chi[0], act[index].bh[1].chi[1],
			 act[index].bh[1].chi[2]);*/
		}
		amp[1] = amp_Start - amp_Step;
	}
	maxT = 0.;
	multi_Match(&params, act, length, dir);
	fprintf(stdout, "%d: %lg\n", maxTindex, maxT);
	maxT = 0.;
	memset(act, 0, length * sizeof(binary_System));
	amp[0] = amp_Start - amp_Step;
	amp[1] = amp_Start - amp_Step;
	lalDebugLevel = 0;
	for (i = 0; i < big_Step; i++) {
		amp[0] += amp_Step;
		for (j = 0; j < big_Step; j++) {
			amp[1] += amp_Step;
			index = i * big_Step + j;
			act[index].F.dec = 0.;
			act[index].F.pol = 0.;
			act[index].F.alpha = 0.;
			act[index].F.gmst = 0.;
			calc_Response_For_Detector(LH, &act[index]);
			act[index].incl = 1.34;
			act[index].dist = 1.;
			act[index].bh[0].m = 3.4;
			act[index].bh[1].m = 3.4;
			convert_Masses(&act[index], FROM_M1M2);
			for (k = 0; k < 2; k++) {
				act[index].bh[k].chi_Amp = amp[k];
				act[index].bh[k].theta = theta[k];
				act[index].bh[k].varphi = varphi[k];
			}
			convert_Spins(&act[index], FROM_THETA_VPHI);
		}
		amp[1] = amp_Start - amp_Step;
	}
	multi_Match(&params, act, length, dir);
	fprintf(stdout, "%d: %lg\n", maxTindex, maxT);
}

void generate(void) {
	char dir[50];
	sprintf(dir, "data_TDK/spinConf/");
	program_Params params;
	memset(&params, 0, sizeof(program_Params));
	params.freq_Sampling = 16380.;
	params.time_Sampling = 1. / params.freq_Sampling;
	params.freq_Initial = 40.;
	params.freq_Final = 0.;
	binary_System min, max;
	init_Binary_System(&min, &max);
	min.bh[0].chi_Amp = max.bh[0].chi_Amp = 0.998;
	min.bh[1].chi_Amp = max.bh[1].chi_Amp = 0.998;
	min.incl = max.incl = 1.34;
	min.dist = max.dist = 1.;
	min.F.dec = max.F.dec = 0.;
	min.F.pol = max.F.pol = 0.;
	min.F.alpha = max.F.alpha = 0.;
	min.F.gmst = max.F.gmst = 0.;
	short i, length = 100;
	binary_System *act = malloc(2 * length * sizeof(binary_System));
	min.bh[0].m = max.bh[0].m = 34.;
	min.bh[1].m = max.bh[1].m = 3.4;
	for (i = 0; i < length; i++) {
		gen_Parameters(&act[i], &min, &max, M1M2, THETA_VPHI);
		convert_Masses(&act[i], FROM_M1M2);
		if (i == 1) {
			printf("%d: % -35.25lg % -35.25lg\n", i, act[i].bh[0].theta,
					act[i].bh[0].varphi);
			printf("%d: % -35.25lg % -35.25lg\n", i, act[i].bh[1].theta,
					act[i].bh[1].varphi);
			fflush(stdout);
			exit(-1);
		}
	}
	multi_Match(&params, act, length, dir);
	for (i = 0; i < length; i++) {
		act[i].bh[0].m = act[i].bh[1].m = 3.4;
		convert_Masses(&act[i], FROM_M1M2);
	}
	multi_Match(&params, act, length, dir);
}
