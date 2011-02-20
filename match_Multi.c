#include "match_Multi.h"

#define M_NUM 4
double calc_Periods(double *per1, double *per2, signalStruct *signal);

void destroySTWave(CoherentGW waveform);

int lalDebugLevel = 0;

double maxT = 0.;
short maxTindex;

extern double theta[4];

extern double varphi[4];

double amp = 0.998;

#define PN1 "twoPN"

void multi_Match(program_Params *params, binary_System *act, long num,
		char dir[50]) {
	assert(params);
	assert(act);
	assert(num>0);
	// saját
	signalStruct sig[M_NUM];
	double *waves[2];
	char file_Name[50];
	FILE *file;
	binary_System sys;
	memset(&sys, 0, sizeof(binary_System));
	// LAL
	static LALStatus status;
	CoherentGW waveform[2];
	SimInspiralTable injParams[2];
	PPNParamStruc ppnParams;
	static RandomInspiralSignalIn randIn;
	memset(&waveform, 0, 2 * sizeof(CoherentGW));
	//	init_Binary_System(&min, &max);
	//	srand(10);
	short j, l;
	long i, k;
	double cshift, sshift, cphi[2], sphi[2], hp;
	double f0, f1;
	for (i = 0; i < num; i++) {
		if ((i + 1) % 10 == 0) {
			fprintf(stdout, "%ld %ld\n", num, i + 1);
		}
		params->index = i + 1;
		//if(i==2)return;
		//		gen_Parameters(&act[i], &min, &max, ETAM, THETA_VPHI);
		memset(&injParams, 0, 2 * sizeof(SimInspiralTable));
		memset(&ppnParams, 0, sizeof(PPNParamStruc));
		injParams[1].mass1 = injParams[0].mass1 = act[i].bh[0].m;
		injParams[1].mass2 = injParams[0].mass2 = act[i].bh[1].m;
		/**//**/
		memcpy(&sys, &act[i], sizeof(binary_System));
		/*
		 sys.bh[0].chi[0] = (sys.bh[1].chi[0] = amp * sin(act[i].incl));
		 sys.bh[0].chi[1] = (sys.bh[1].chi[1] = 0.);
		 sys.bh[0].chi[2] = (sys.bh[1].chi[2] = amp * cos(act[i].incl));
		 *//*
		 sys.bh[0].chi[0] = -(sys.bh[1].chi[0] = -amp * cos(act[i].incl));
		 sys.bh[0].chi[1] = -(sys.bh[1].chi[1] = 0.);
		 sys.bh[0].chi[2] = -(sys.bh[1].chi[2] = amp * sin(act[i].incl));
		 *//*
		sys.bh[0].chi[0] = -(sys.bh[1].chi[0] = amp * cos(act[i].incl));
		sys.bh[0].chi[1] = -(sys.bh[1].chi[1] = 0.);
		sys.bh[0].chi[2] = -(sys.bh[1].chi[2] = amp * sin(act[i].incl));
		*//*
		 sys.bh[0].chi[0] = -(sys.bh[1].chi[0] = -amp * cos(act[i].incl));
		 sys.bh[0].chi[1] = -(sys.bh[1].chi[1] = 0.);
		 sys.bh[0].chi[2] = -(sys.bh[1].chi[2] = amp * sin(act[i].incl));
		 *//*
		 sys.bh[0].chi[0] = amp * cos(act[i].incl);
		 sys.bh[0].chi[1] = 0.;
		 sys.bh[0].chi[2] = -amp * sin(act[i].incl);
		 sys.bh[1].chi[0] = 0.;
		 sys.bh[1].chi[1] = amp;
		 sys.bh[1].chi[2] = 0.;
		 convert_Spins(&sys, FROM_XYZ);
		 *//*
		 sys.bh[0].chi_Amp = sys.bh[1].chi_Amp = amp;
		 sys.bh[0].kappa = M_PI / 6.;
		 sys.bh[0].psi = 0.;
		 sys.bh[1].kappa = M_PI / 6.;
		 sys.bh[1].psi = M_PI;
		 convert_Spins(&sys, FROM_KAPPA_PSI);*/
		/*
		injParams[0].spin1x = sys.bh[0].chi[0];
		injParams[0].spin1y = sys.bh[0].chi[1];
		injParams[0].spin1z = sys.bh[0].chi[2];
		injParams[0].spin2x = sys.bh[1].chi[0];
		injParams[0].spin2y = sys.bh[1].chi[1];
		injParams[0].spin2z = sys.bh[1].chi[2];
		*/
		 injParams[0].spin1x = act[i].bh[0].chi[0];
		 injParams[0].spin1y = act[i].bh[0].chi[1];
		 injParams[0].spin1z = act[i].bh[0].chi[2];
		 injParams[0].spin2x = act[i].bh[1].chi[0];
		 injParams[0].spin2y = act[i].bh[1].chi[1];
		 injParams[0].spin2z = act[i].bh[1].chi[2];
		 /**/
		injParams[1].spin1x = act[i].bh[0].chi[0];
		injParams[1].spin1y = act[i].bh[0].chi[1];
		injParams[1].spin1z = act[i].bh[0].chi[2];
		injParams[1].spin2x = act[i].bh[1].chi[0];
		injParams[1].spin2y = act[i].bh[1].chi[1];
		injParams[1].spin2z = act[i].bh[1].chi[2];
		injParams[1].inclination = injParams[0].inclination = act[i].incl;
		double freq_Min = injParams[0].f_lower = 40.;
		injParams[1].f_final = injParams[0].f_final = 0.;
		injParams[1].distance = injParams[0].distance = act[i].dist;
		injParams[1].polarization = injParams[0].polarization = 0;
		injParams[1].coa_phase = injParams[0].coa_phase = act[i].coaPhase = 0.;
		ppnParams.deltaT = 1. / params->freq_Sampling;
		injParams[1].f_lower = injParams[0].f_lower = params->freq_Initial;
		snprintf(injParams[0].waveform,
				LIGOMETA_WAVEFORM_MAX * sizeof(CHAR),
				"SpinQuadTaylor"PN1"SOSSQM");
		snprintf(injParams[1].waveform,
				LIGOMETA_WAVEFORM_MAX * sizeof(CHAR), "SpinTaylor"PN1"SS");
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
		long min_Length = waveform[shorter].f->data->length;
		long max_Length = waveform[!shorter].f->data->length;
		act[i].coaTime = (waveform[0].f->data->length - 1)
				* params->time_Sampling;
		f0 = waveform[0].f->data->data[waveform[0].f->data->length - 1];
		f1 = waveform[1].f->data->data[waveform[1].f->data->length - 1];
		create_Signal_Struct(&sig[0], waveform[!shorter].f->data->length);
		create_Signal_Struct(&sig[1], waveform[!shorter].f->data->length);
		create_Signal_Struct(&sig[2], waveform[!shorter].f->data->length);
		create_Signal_Struct(&sig[3], waveform[!shorter].f->data->length);
		waves[0] = malloc(max_Length * sizeof(double));
		waves[1] = malloc(max_Length * sizeof(double));
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
				waves[j][k] = act[i].F.F[0] * (waveform[j].a->data->data[2 * k]
						* cshift * cphi[0] - waveform[j].a->data->data[2 * k
						+ 1] * sshift * sphi[0]) + act[i].F.F[1]
						* (waveform[j].a->data->data[2 * k] * sshift * cphi[0]
								+ waveform[j].a->data->data[2 * k + 1] * cshift
										* sphi[0]);
				for (l = 0; l < M_NUM; l++) {
					sig[l].signal[2 * j][k] = act[i].F.F[0] * //
							(waveform[j].a->data->data[2 * k] * cshift * //
									cphi[0] - waveform[j].a->data->data[2 * k
									+ 1] * sshift * sphi[0]) + act[i].F.F[1] * //
							(waveform[j].a->data->data[2 * k] * sshift * //
									cphi[0] + waveform[j].a->data->data[2 * k
									+ 1] * cshift * sphi[0]);
					sig[l].signal[2 * j + 1][k] = act[i].F.F[0] * //
							(waveform[j].a->data->data[2 * k] * cshift * //
									cphi[1] - waveform[j].a->data->data[2 * k
									+ 1] * sshift * sphi[1]) + act[i].F.F[1] * //
							(waveform[j].a->data->data[2 * k] * sshift * //
									cphi[1] + waveform[j].a->data->data[2 * k
									+ 1] * cshift * sphi[1]);
				}
			}
		}
		params->periodsD = calc_Periods(&params->periods[0],
				&params->periods[1], &sig[0]);
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
		params->match_Typ = match_Typical(&sig[0], minfr, maxfr);
		//printf("%ld %d: %lg %lg\n", i, maxTindex, maxT, params->match_TypT);
		if (params->match_Typ > maxT && params->match_Typ < 1.) {
			maxT = params->match_Typ;
			maxTindex = i;
			printf("%d: % 20.15lg %lg\n", maxTindex, maxT, act[i].bh[0].chi_Amp);
		}
		calc_Overlap_Time(&params->match_Best, &params->match_Worst, &sig[3],
				minfr, maxfr);
		sprintf(file_Name, "%s%dgen%04ld.dat", dir, (int) act[i].bh[0].m, i);
		file = fopen(file_Name, "w");
		print_Binary_System(&act[i], params, file, f0, f1);
		fprintf(file, "#%s%16s%16s\n", "time", "hSQT", "hST");
		for (k = 0; k < min_Length; k++) {
			fprintf(file, PREC" X "PREC" X "PREC"\n", k * ppnParams.deltaT,
					waves[0][k], waves[1][k]);
		}
		for (; k < max_Length; k++) {
			if (shorter) {
				fprintf(file, PREC" X "PREC"\n", k * ppnParams.deltaT,
						waves[0][k]);
			} else {
				fprintf(file, PREC" X ""%-14.8s "" X "PREC"\n", k
						* ppnParams.deltaT, "", waves[1][k]);
			}
		}
		fclose(file);
		XLALSQTPNDestroyCoherentGW(&waveform[0]);
		//XLALSQTPNDestroyCoherentGW(&waveform[1]);
		destroySTWave(waveform[1]);
		destroy_Signal_Struct(&sig[0]);
		destroy_Signal_Struct(&sig[1]);
		destroy_Signal_Struct(&sig[2]);
		destroy_Signal_Struct(&sig[3]);
		free(waves[0]);
		free(waves[1]);
		XLALFree(randIn.psd.data);
	}
	hp = params->freq_Sampling;
	memset(params, 0, sizeof(program_Params));
	params->freq_Sampling = hp;
	params->time_Sampling = 1. / params->freq_Sampling;
	params->freq_Initial = 40.;
	//		destroySTWave(waveform[0]);
	//		destroySTWave(waveform[1]);
	//destroySTWave(waveform[1]);
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

double calc_Periods(double *per1, double *per2, signalStruct *signal) {
	assert(signal);
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

binary_System* fill_TDK(long *length, binary_System *minta, cover cov) {
	binary_System *sys;
	if (cov == SPECIAL) {
		*length = 6;
		double amp = 0.998;
		sys = malloc(*length * sizeof(binary_System));
		short i;
		for (i = 0; i < *length; i++) {
			memcpy(&sys[i], minta, sizeof(binary_System));
		}
		i = 0;
		// pályasíkra merőleges konfiguráció
		sys[i].bh[0].chi[0] = sys[i].bh[1].chi[0] = amp * sin(sys->incl);
		sys[i].bh[0].chi[1] = sys[i].bh[1].chi[1] = 0.;
		sys[i].bh[0].chi[2] = sys[i].bh[1].chi[2] = amp * cos(sys->incl);
		convert_Spins(&sys[i], FROM_XYZ);
		i++;
		// pályasíkos konfigurációk
		sys[i].bh[0].chi[0] = -(sys[i].bh[1].chi[0] = -amp * cos(sys->incl));
		sys[i].bh[0].chi[1] = -(sys[i].bh[1].chi[1] = 0.);
		sys[i].bh[0].chi[2] = -(sys[i].bh[1].chi[2] = amp * sin(sys->incl));
		convert_Spins(&sys[i], FROM_XYZ);
		i++;
		// pályasíkos konfigurációk2
		sys[i].bh[0].chi[0] = -(sys[i].bh[1].chi[0] = amp * cos(sys->incl));
		sys[i].bh[0].chi[1] = -(sys[i].bh[1].chi[1] = 0.);
		sys[i].bh[0].chi[2] = -(sys[i].bh[1].chi[2] = amp * sin(sys->incl));
		convert_Spins(&sys[i], FROM_XYZ);
		i++;
		// merőleges konfiguráció
		sys[i].bh[0].chi[0] = amp * cos(sys->incl);
		sys[i].bh[0].chi[1] = 0.;
		sys[i].bh[0].chi[2] = -amp * sin(sys->incl);
		sys[i].bh[1].chi[0] = 0.;
		sys[i].bh[1].chi[1] = amp;
		sys[i].bh[1].chi[2] = 0.;
		convert_Spins(&sys[i], FROM_XYZ);
		i++;
		// tengely-szimmetrikus konfiguráció
		sys[i].bh[0].chi_Amp = sys[i].bh[1].chi_Amp = amp;
		sys[i].bh[0].kappa = M_PI / 6.;
		sys[i].bh[0].psi = 0.;
		sys[i].bh[1].kappa = M_PI / 6.;
		sys[i].bh[1].psi = M_PI;
		convert_Spins(&sys[i], FROM_KAPPA_PSI);
		i++;
		// pont-szimmetrikus konfiguráció
		sys[i].bh[0].chi_Amp = sys[i].bh[1].chi_Amp = amp;
		sys[i].bh[0].kappa = M_PI / 6.;
		sys[i].bh[0].psi = 0.;
		sys[i].bh[1].kappa = M_PI - M_PI / 6.;
		sys[i].bh[1].psi = M_PI;
		convert_Spins(&sys[i], FROM_KAPPA_PSI);
	} else {
		sys = malloc(*length * sizeof(binary_System));
		double step = 0.001;
		double min = 0.998;
		double max = 1.200;
		*length = (long) ((max - min) / step);
		printf("%ld\n", *length);
		double amp = min;
		sys = malloc(*length * sizeof(binary_System));
		short i;
		for (i = 0; i < *length; i++) {
			memcpy(&sys[i], minta, sizeof(binary_System));
			/*
			 sys[i].bh[0].chi[0] = sys[i].bh[1].chi[0] = amp * sin(sys->incl);
			 sys[i].bh[0].chi[1] = sys[i].bh[1].chi[1] = 0.;
			 sys[i].bh[0].chi[2] = sys[i].bh[1].chi[2] = amp * cos(sys->incl);
			 *//*
			 sys[i].bh[0].chi[0]
			 = -(sys[i].bh[1].chi[0] = -amp * cos(sys->incl));
			 sys[i].bh[0].chi[1] = -(sys[i].bh[1].chi[1] = 0.);
			 sys[i].bh[0].chi[2] = -(sys[i].bh[1].chi[2] = amp * sin(sys->incl));
			 */
			sys[i].bh[0].chi[0] = -(sys[i].bh[1].chi[0] = amp * cos(sys->incl));
			sys[i].bh[0].chi[1] = -(sys[i].bh[1].chi[1] = 0.);
			sys[i].bh[0].chi[2] = -(sys[i].bh[1].chi[2] = amp * sin(sys->incl));
			/*
			 sys[i].bh[0].chi[0] = amp * cos(sys->incl);
			 sys[i].bh[0].chi[1] = 0.;
			 sys[i].bh[0].chi[2] = -amp * sin(sys->incl);
			 sys[i].bh[1].chi[0] = 0.;
			 sys[i].bh[1].chi[1] = amp;
			 sys[i].bh[1].chi[2] = 0.;
			 */
			convert_Spins(&sys[i], FROM_XYZ);
			/*sys[i].bh[0].chi_Amp = sys[i].bh[1].chi_Amp = amp;
			 sys[i].bh[0].kappa = M_PI / 6.;
			 sys[i].bh[0].psi = 0.;
			 sys[i].bh[1].kappa = M_PI / 6.;
			 sys[i].bh[1].psi = M_PI;
			 convert_Spins(&sys[i], FROM_KAPPA_PSI);
			 */amp += step;
		}
	}
	return sys;
}
