#include "match_Multi.h"

#define M_NUM 4
double calc_Periods(double *per1, double *per2, signalStruct *signal);

void destroySTWave(CoherentGW waveform);

void multi_Match(program_Params *params, binary_System *act, long num, char dir[50]) {
	// saját
	signalStruct sig[M_NUM];
	char file_Name[2][50];
	FILE *file[2];
	// LAL
	static LALStatus status;
	CoherentGW waveform[2];
	SimInspiralTable injParams[2];
	PPNParamStruc ppnParams;
	static RandomInspiralSignalIn randIn;
	memset(&waveform, 0, 2 * sizeof(CoherentGW));
	memset(&injParams, 0, 2 * sizeof(SimInspiralTable));
	memset(&ppnParams, 0, sizeof(PPNParamStruc));
//	init_Binary_System(&min, &max);
//	srand(10);
	short j,l;
	long i,k;
	double cshift, sshift, cphi[2], sphi[2], hp, hc;
	for (i = 0; i < num; i++) {
		printf("%ld %ld\n", num, i + 1);
		params->index = i + 1;
//		gen_Parameters(&act[i], &min, &max, ETAM, THETA_VPHI);
		injParams[0].mass1 = act[i].bh[0].m;
		injParams[0].mass2 = act[i].bh[1].m;
		injParams[0].spin1x = act[i].bh[0].chi[0];// = 0.;//998 * sin(1.43);// = 0.0790765695;
		injParams[0].spin1y = act[i].bh[0].chi[1];// = 0.998;// = 0.1589679868;
		injParams[0].spin1z = act[i].bh[0].chi[2];// = 0.;//998 * cos(1.43);// = 0.9820794649;
		injParams[0].spin2x = act[i].bh[1].chi[0];// = 0.;//998 * sin(1.43);// = 0.118388124;
		injParams[0].spin2y = act[i].bh[1].chi[1];// = 0.998;// = -0.1715059157;
		injParams[0].spin2z = act[i].bh[1].chi[2];// = 0.;//998 * cos(1.43);// = 0.9759989616;
		injParams[0].qmParameter1 = 1.;
		injParams[0].qmParameter2 = 1.;
		injParams[0].inclination = act[i].incl;
		double freq_Min = injParams[0].f_lower = 40.;
		injParams[0].f_final = 0.;
		injParams[0].distance = act[i].dist;
		injParams[0].polarization = 0;
		injParams[0].coa_phase = act[i].coaPhase;
		ppnParams.deltaT = 1./params->freq_Sampling;
		injParams[0].f_lower = params->freq_Initial;
		memcpy(&injParams[1], &injParams[0], sizeof(SimInspiralTable));
		LALSnprintf(injParams[0].waveform,
				LIGOMETA_WAVEFORM_MAX * sizeof(CHAR), "SpinQuadTaylor"PN"SS");
		LALSnprintf(injParams[1].waveform,
				LIGOMETA_WAVEFORM_MAX * sizeof(CHAR), "SpinQuadTaylor"PN"ALL");
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
		long shorter = waveform[0].f->data->length
				< waveform[1].f->data->length ? 0 : 1;
		long length = waveform[!shorter].f->data->length;
		waveform[!shorter].f->data->length = length;
		create_Signal_Struct(&sig[0], waveform[!shorter].f->data->length);
		create_Signal_Struct(&sig[1], waveform[!shorter].f->data->length);
		create_Signal_Struct(&sig[2], waveform[!shorter].f->data->length);
		create_Signal_Struct(&sig[3], waveform[!shorter].f->data->length);
		randIn.psd.length = waveform[!shorter].f->data->length;
		double df = 1. / ppnParams.deltaT / randIn.psd.length;
		randIn.psd.data = (REAL8*) LALMalloc(sizeof(REAL8) * randIn.psd.length);
		LALNoiseSpectralDensity(&status, &randIn.psd, &LALLIGOIPsd, df);
		for (j = 0; j < randIn.psd.length; j++) {
			sig[0].psd[j] = sig[1].psd[j] = sig[2].psd[j] = sig[3].psd[j]
					= randIn.psd.data[j];
		}
		for (j = 0; j < 2; j++) {
			for (k = 0; k < waveform[j].f->data->length; k++) {
				cshift = cos(waveform[j].shift->data->data[k]);
				sshift = sin(waveform[j].shift->data->data[k]);
				cphi[0] = cos(waveform[j].phi->data->data[k]);
				sphi[0] = sin(waveform[j].phi->data->data[k]);
				cphi[1] = cos(waveform[j].phi->data->data[k] + M_PI / 2.0);
				sphi[1] = sin(waveform[j].phi->data->data[k] + M_PI / 2.0);
/*				hp = waveform[j].a->data->data[2 * k] * cshift * cphi
						- waveform[j].a->data->data[2 * k + 1] * sshift * sphi;
				hc = waveform[j].a->data->data[2 * k] * sshift * cphi
						+ waveform[j].a->data->data[2 * k + 1] * cshift * sphi;*/
				for (l = 0; l < M_NUM; l++) {
					sig[l].signal[2 * j][k] = act[i].F.F[0] * (
						waveform[j].a->data->data[2 * k]     * cshift * cphi[0] -
						waveform[j].a->data->data[2 * k + 1] * sshift * sphi[0]) + act[i].F.F[1] * (
						waveform[j].a->data->data[2 * k]     * sshift * cphi[0] +
						waveform[j].a->data->data[2 * k + 1] * cshift * sphi[0]);
					sig[l].signal[2 * j + 1][k] = act[i].F.F[0] * (
						waveform[j].a->data->data[2 * k]     * cshift * cphi[1] -
						waveform[j].a->data->data[2 * k + 1] * sshift * sphi[1]) + act[i].F.F[1] * (
						waveform[j].a->data->data[2 * k]     * sshift * cphi[1] +
						waveform[j].a->data->data[2 * k + 1] * cshift * sphi[1]);
				}
			}
		}
		params->periodsD = calc_Periods(&params->periods[0], &params->periods[1], &sig[0]);
		act[i].coaTime[0] = (waveform[0].f->data->length - 1) * params->time_Sampling;
		act[i].coaTime[1] = (waveform[1].f->data->length - 1) * params->time_Sampling;
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
		//		printf("Simp = "PREC"\n", match_Simple(&simp, minfr, maxfr));fflush(stdout);
		params->match_Typ = match_Typical(&sig[0], minfr, maxfr, NONE);
		params->match_TypT = match_Typical(&sig[1], minfr, maxfr, TIME);
		calc_Overlap(&params->match_Best, &params->match_Worst, &sig[2], minfr, maxfr);
		calc_Overlap_Time(&params->match_BestT, &params->match_WorstT, &sig[3], minfr, maxfr);
//		print_Binary_System(&act[i], params, stdout, act[i].coaTime[0]);
		for (j = 0; j < 2; j++) {
			sprintf(file_Name[j], "%sgen%04ld_%d.dat", dir, i, j);
			file[j] = fopen(file_Name[j], "w");
			print_Binary_System(&act[i], params, file[j], act[i].coaTime[j]);
			fprintf(file[j], "#%s%16s%16s%16s\n", "time", "h", "h_+", "h_x");
			for (k = 0; k < waveform[j].f->data->length; k++) {
				fprintf(file[j], PREC PREC PREC PREC"\n", k * ppnParams.deltaT,
						sig[0].signal[2 * j][k] * act[i].F.F[0] + sig[0].signal[2
								* j + 1][k] * act[i].F.F[1],
						sig[0].signal[2 * j][k], sig[0].signal[2 * j + 1][k]);
			}
			fclose(file[j]);
		}
		hp = params->freq_Sampling;
		memset(params, 0, sizeof(program_Params));
		params->freq_Sampling = hp;
		params->time_Sampling = 1. / params->freq_Sampling;
		params->freq_Initial = 40.;
		XLALSQTPNDestroyCoherentGW(&waveform[0]);
		XLALSQTPNDestroyCoherentGW(&waveform[1]);
//		destroySTWave(waveform[0]);
//		destroySTWave(waveform[1]);
		destroy_Signal_Struct(&sig[0]);
		destroy_Signal_Struct(&sig[1]);
		destroy_Signal_Struct(&sig[2]);
		destroy_Signal_Struct(&sig[3]);
		//destroySTWave(waveform[1]);
		XLALFree(randIn.psd.data);
	}
}

void destroySTWave(CoherentGW waveform) {
	if (waveform.f->data)
		XLALDestroyREAL4Vector(waveform.f->data);
	if (waveform.f)
		LALFree(waveform.f);
//	if (waveform.shift->data)
//		XLALDestroyREAL4Vector(waveform.shift->data);
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

double calc_Periods(double *per1, double *per2, signalStruct *signal) {
	double prev, act;
	*per1 = *per2 = 0;
	long i;
	prev = signal->signal[H1][0] * signal->signal[H1][1];
	for (i = 1; i < signal->size; i++) {
		act = signal->signal[H1][i] * signal->signal[H1][i + 1];
		if (act <= 0. && prev > 0.) {
			(*per1)++;
		}
		prev = act;
	}
	prev = signal->signal[H2][0] * signal->signal[H2][1];
	for (i = 1; i < signal->size; i++) {
		act = signal->signal[H2][i] * signal->signal[H2][i + 1];
		if (act <= 0. && prev > 0.) {
			(*per2)++;
		}
		prev = act;
	}
	(*per1) /= 2.;
	(*per2) /= 2.;
	return (*per1 - *per2) / (*per1);
}

