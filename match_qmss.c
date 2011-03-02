/**
 * @file match_qmss.c
 *
 * @date Feb 25, 2011
 * @author vereb
 */

#include "match_qmss.h"

short is_First;///<a

void read_Program_Parameters(Program_Parameters *parameters, System_Parameters *params,
		char *file_Name) {
	assert(parameters);
	assert(params);
	assert(file_Name);
	FILE *file = fopen(file_Name, "r");
	fscanf(file, "%ld\n", &parameters->number_Of_Runs);
	fscanf(file, "%hd\n", &parameters->precision);
	parameters->width_Of_Number = parameters->precision + EXTRA_CHARACTERS;
	fscanf(file, "%hd\n", &parameters->precision_To_Plot);
	parameters->width_Of_Number_To_Plot = parameters->precision_To_Plot + EXTRA_CHARACTERS;
	fscanf(file, "%s\n", parameters->folder);
	fscanf(file, "%lg\n", &params->min_Match);
	fscanf(file, "%lg\n", &params->max_Spin);
	fscanf(file, "%lg\n", &params->spin_Step);
	fscanf(file, "%lg\n", &params->freq_Sampling);
	params->time_Sampling = 1. / params->freq_Sampling;
	fscanf(file, "%lg\n", &params->freq_Initial);
	fscanf(file, "%lg\n", &params->delta_Length);
	fscanf(file, "%s\n", params->approx[0]);
	fscanf(file, "%s\n", params->phase[0]);
	fscanf(file, "%hd\n", &params->amp_Code);
	fscanf(file, "%s\n", params->spin[0]);
	fscanf(file, "%s\n", params->approx[1]);
	fscanf(file, "%s\n", params->phase[1]);
	fscanf(file, "%hd\n", &params->amp_Code);
	fscanf(file, "%s\n", params->spin[1]);
	fclose(file);
}

void read_Parameters(binary_System *parameters, char *file_Name) {
	assert(parameters);
	assert(file_Name);
	FILE *file = fopen(file_Name, "r");
	fscanf(file, "%lg %lg\n", &parameters[0].M, &parameters[1].M);
	fscanf(file, "%lg %lg\n", &parameters[0].eta, &parameters[1].eta);
	fscanf(file, "%lg %lg\n", &parameters[0].bh[0].m, &parameters[1].bh[0].m);
	fscanf(file, "%lg %lg\n", &parameters[0].bh[1].m, &parameters[1].bh[1].m);
	fscanf(file, "%lg %lg\n", &parameters[0].bh[0].chi_Amp, &parameters[1].bh[0].chi_Amp);
	fscanf(file, "%lg %lg\n", &parameters[0].bh[1].chi_Amp, &parameters[1].bh[1].chi_Amp);
	fscanf(file, "%lg %lg\n", &parameters[0].bh[0].kappa, &parameters[1].bh[0].kappa);
	parameters[0].bh[0].kappa *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
	parameters[1].bh[0].kappa *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
	fscanf(file, "%lg %lg\n", &parameters[0].bh[1].kappa, &parameters[1].bh[1].kappa);
	parameters[0].bh[1].kappa *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
	parameters[1].bh[1].kappa *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
	fscanf(file, "%lg %lg\n", &parameters[0].bh[0].psi, &parameters[1].bh[0].psi);
	parameters[0].bh[0].psi *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
	parameters[1].bh[0].psi *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
	fscanf(file, "%lg %lg\n", &parameters[0].bh[1].psi, &parameters[1].bh[1].psi);
	parameters[0].bh[1].psi *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
	parameters[1].bh[1].psi *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
	fscanf(file, "%lg %lg\n", &parameters[0].incl, &parameters[1].incl);
	fscanf(file, "%lg %lg\n", &parameters[0].dist, &parameters[1].dist);
	fscanf(file, "%lg %lg\n", &parameters[0].F.pol, &parameters[1].F.pol);
	parameters[0].F.pol *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
	parameters[1].F.pol *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
	fscanf(file, "%lg %lg\n", &parameters[0].F.alpha, &parameters[1].F.alpha);
	fscanf(file, "%lg %lg\n", &parameters[0].F.dec, &parameters[1].F.dec);
	parameters[0].F.dec *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
	parameters[1].F.dec *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
	fscanf(file, "%lg %lg\n", &parameters[0].F.gmst, &parameters[1].F.gmst);
	fclose(file);
}

void run_Algorithm(Program_Parameters *program_Parameters, System_Parameters *parameters,
		binary_System *limits) {
	assert(program_Parameters);
	assert(parameters);
	assert(limits);
	assert(program_Parameters->number_Of_Runs > 0);
	srand(86);
	short is_Good;
	char temp[FILE_NAME_LENGTH];
	sprintf(temp, "%s", program_Parameters->folder);
	for (long i = 0; i < program_Parameters->number_Of_Runs;) {
		generate_Parameters(parameters, limits);
		is_Good = incrementing_Spins(program_Parameters, parameters);
		if (is_Good) {
			i++;
		}
		sprintf(program_Parameters->folder, "%s", temp);
	}
}

void generate_Parameters(System_Parameters *parameters, binary_System *limits) {
	assert(parameters);
	assert(limits);
	gen_Parameters(&parameters->system[0], &limits[0], &limits[1], ETAM, KAPPA_PSI);
	memcpy(&parameters->system[1], &parameters->system[0], sizeof(binary_System));
}

short incrementing_Spins(Program_Parameters *prog, System_Parameters* parameters) {
	assert(prog);
	assert(parameters);
	char temp[FILE_NAME_LENGTH];
	sprintf(temp, "%s", prog->folder);
	is_First = 1;
	short is_Good;
	parameters->critical_Match = 0.0;
	for (; parameters->system[MOD_SPIN_INDEX].bh[0].chi_Amp < parameters->max_Spin; parameters->system[MOD_SPIN_INDEX].bh[0].chi_Amp
			+= parameters->spin_Step) {
		for (; parameters->system[MOD_SPIN_INDEX].bh[1].chi_Amp < parameters->max_Spin; parameters->system[MOD_SPIN_INDEX].bh[1].chi_Amp
				+= parameters->spin_Step) {
			is_Good = calc_Matches_For_ParameterPair(prog, parameters);
			if (is_First) {
				if (!is_Good) {
					return is_Good;
				}
				parameters->critical_Match = parameters->match_Minimax;
				is_First = 0;
			}
			sprintf(prog->folder, "%s", temp);
		}
	}
	return is_Good;
}

inline void increment_Spin_Of_Binary_System(binary_System *system, double step) {
	assert(system);
	assert(step>0.);
	system->bh[0].chi_Amp += step;
	system->bh[1].chi_Amp += step;
}

inline void increment_Spins(System_Parameters* parameters) {
	assert(parameters);
	increment_Spin_Of_Binary_System(&parameters->system[0], parameters->spin_Step);
	increment_Spin_Of_Binary_System(&parameters->system[1], parameters->spin_Step);
}

short calc_Matches_For_ParameterPair(Program_Parameters *prog, System_Parameters *parameters) {
	assert(prog);
	assert(parameters);
	static LALParameters lalparams;
	//	double f0, f1;
	signalStruct sig;
	initLALParameters(&lalparams, parameters);
	for (short i = 0; i < 2; i++) {
		memset(&lalparams.status, 0, sizeof(LALStatus));
		LALGenerateInspiral(&lalparams.status, &lalparams.waveform[i], &lalparams.injParams[i],
				&lalparams.ppnParams);
		if (lalparams.status.statusCode) {
			fprintf(stderr, "%d: LALSQTPNWaveformTest: error generating waveform\n", i);
			XLALSQTPNDestroyCoherentGW(&lalparams.waveform[0]);
			return 1;
		}
		parameters->system[i].coaPhase
				= lalparams.waveform[i].phi->data->data[lalparams.waveform[i].phi->data->length - 1];
		parameters->system[i].coaTime = (lalparams.waveform[i].f->data->length - 1)
				* parameters->time_Sampling;
	}
	parameters->shorter = lalparams.shorter = lalparams.waveform[0].f->data->length
			< lalparams.waveform[1].f->data->length ? 0 : 1;
	parameters->min_Length = lalparams.min_Length
			= lalparams.waveform[lalparams.shorter].f->data->length;
	parameters->max_Length = lalparams.max_Length
			= lalparams.waveform[!lalparams.shorter].f->data->length;
	parameters->freq_Max = (lalparams.injParams[0].f_final + lalparams.injParams[1].f_final) / 2.;
	parameters->freq_Step = 1. / (lalparams.ppnParams.deltaT * lalparams.max_Length);
	create_Signal_Struct(&sig, lalparams.waveform[!lalparams.shorter].f->data->length);
	createPSD(&lalparams, &sig);
	for (short i = 0; i < 2; i++) {
		for (long j = 0; j < lalparams.waveform[i].f->data->length; j++) {
			sig.signal[2 * i][j] = lalparams.waveform[i].h->data->data[2 * j];
			sig.signal[2 * i + 1][j] = lalparams.waveform[i].h->data->data[2 * j + 1];
		}
	}
	double fr = 0.;
	long minfr = 0, maxfr = 0;
	while (fr < parameters->freq_Min) {
		fr += parameters->freq_Step;
		maxfr = ++minfr;
	}
	while (fr < parameters->freq_Max) {
		fr += parameters->freq_Step;
		maxfr++;
	}
	parameters->match_Typ = parameters->match_Best = parameters->match_Minimax = 0.0;
	calc_Matches(&sig, minfr, maxfr, &parameters->match_Typ, &parameters->match_Best,
			&parameters->match_Minimax);
	printf("% 11.6lg% 11.6lg% 11.6lg\n", parameters->match_Typ, parameters->match_Best,
			parameters->match_Minimax);
	if (parameters->match_Minimax > parameters->min_Match) {
		write_Waves_To_Files(prog, parameters, &sig);
	}
	XLALSQTPNDestroyCoherentGW(&lalparams.waveform[0]);
	XLALSQTPNDestroyCoherentGW(&lalparams.waveform[1]);
	destroy_Signal_Struct(&sig);
	XLALFree(lalparams.randIn.psd.data);
	if (parameters->match_Minimax > parameters->min_Match) {
		return 1;
	}
	return 0;
}

void initLALParameters(LALParameters *lalparams, System_Parameters *parameters) {
	assert(lalparams);
	assert(parameters);
	memset(&lalparams->waveform, 0, 2 * sizeof(CoherentGW));
	memset(&lalparams->injParams, 0, 2 * sizeof(SimInspiralTable));
	memset(&lalparams->ppnParams, 0, sizeof(PPNParamStruc));
	memset(&lalparams->waveform, 0, 2 * sizeof(CoherentGW));
	lalparams->ppnParams.deltaT = 1. / parameters->freq_Sampling;
	parameters->freq_Min = 40.;
	for (short i = 0; i < 2; i++) {
		lalparams->injParams[i].mass1 = parameters->system[i].bh[0].m;
		lalparams->injParams[i].mass2 = parameters->system[i].bh[1].m;
		lalparams->injParams[i].spin1x = parameters->system[i].bh[0].chi[0];
		lalparams->injParams[i].spin1y = parameters->system[i].bh[0].chi[1];
		lalparams->injParams[i].spin1z = parameters->system[i].bh[0].chi[2];
		lalparams->injParams[i].spin2x = parameters->system[i].bh[1].chi[0];
		lalparams->injParams[i].spin2y = parameters->system[i].bh[1].chi[1];
		lalparams->injParams[i].spin2z = parameters->system[i].bh[1].chi[2];
		lalparams->injParams[i].inclination = parameters->system[i].incl;
		lalparams->injParams[i].f_lower = parameters->freq_Min;
		lalparams->injParams[i].distance = parameters->system[i].dist;
		lalparams->injParams[i].coa_phase = parameters->system[i].coaPhase = 0.;
		lalparams->injParams[i].f_lower = parameters->freq_Initial;
		lalparams->ppnParams.deltaT = 1. / parameters->freq_Sampling;
		lalparams->injParams[i].amp_order = parameters->amp_Code;
		snprintf(lalparams->injParams[i].waveform, LIGOMETA_WAVEFORM_MAX * sizeof(CHAR),
				"%s%s%s", parameters->approx[i], parameters->phase[i],
				parameters->spin[i]);
	}
}

void createPSD(LALParameters *lalparams, signalStruct *sig) {
	assert(lalparams);
	assert(sig);
	lalparams->randIn.psd.length = lalparams->max_Length;
	double df = 1. / lalparams->ppnParams.deltaT / lalparams->randIn.psd.length;
	lalparams->randIn.psd.data = (REAL8*) LALMalloc(sizeof(REAL8) * lalparams->randIn.psd.length);
	LALNoiseSpectralDensity(&lalparams->status, &lalparams->randIn.psd, &LALLIGOIPsd, df);
	for (long j = 0; j < lalparams->randIn.psd.length; j++) {
		sig->psd[j] = lalparams->randIn.psd.data[j];
	}
}

void write_Waves_To_Files(Program_Parameters *prog, System_Parameters *parameters,
		signalStruct *sig) {
	assert(prog);
	assert(parameters);
	assert(sig);
	static short index[3] = { 0, 0, 0 };
	char temp[FILE_NAME_LENGTH];
	sprintf(temp, "%s", prog->folder);
	if (parameters->match_Minimax > parameters->critical_Match && (parameters->max_Length
			- parameters->min_Length) * parameters->time_Sampling < parameters->delta_Length) {
		sprintf(prog->folder, "%s/best", temp);
		write_Wave_To_File(prog, parameters, sig, index[0]);
		index[0]++;
	} else if (parameters->match_Minimax > parameters->critical_Match) {
		sprintf(prog->folder, "%s/match", temp);
		write_Wave_To_File(prog, parameters, sig, index[1]);
		index[1]++;
	}
}

void write_Wave_To_File(Program_Parameters *prog, System_Parameters *parameters, signalStruct *sig,
		short index) {
	assert(prog);
	assert(parameters);
	assert(sig);
	char file_Name[FILE_NAME_LENGTH];
	static char temp[FILE_NAME_LENGTH];
	static char text[FILE_NAME_LENGTH];
	sprintf(temp, "%%%d.%dlg\t", prog->width_Of_Number, prog->precision);
	sprintf(text, "%s%s%s", temp, temp, temp);
	sprintf(file_Name, "%s/data.txt", prog->folder);
	FILE *file = sfopen(file_Name, "a");
	if (file) {
		fprintf(file, text, parameters->system[0].M, parameters->system[0].eta,
				parameters->system[0].chirpM);
		fprintf(file, text, parameters->system[0].mu, parameters->system[0].bh[0].m,
				parameters->system[0].bh[1].m);
		fprintf(file, text, parameters->system[0].incl, parameters->system[0].dist,
				parameters->system[0].F.pol);
		fprintf(file, text, parameters->system[0].F.alpha, parameters->system[0].F.dec,
				parameters->system[0].F.gmst);
		fprintf(file, text, parameters->system[0].bh[0].chi_Amp, parameters->system[0].bh[0].kappa,
				parameters->system[0].bh[0].psi);
		fprintf(file, text, parameters->system[0].bh[1].chi_Amp, parameters->system[0].bh[1].kappa,
				parameters->system[0].bh[1].psi);
		fprintf(file, text, parameters->system[0].bh[0].chi[0], parameters->system[0].bh[0].chi[1],
				parameters->system[0].bh[0].chi[2]);
		fprintf(file, text, parameters->system[0].bh[1].chi[0], parameters->system[0].bh[1].chi[1],
				parameters->system[0].bh[1].chi[2]);
		fprintf(file, text, parameters->system[0].bh[0].theta, parameters->system[0].bh[0].varphi,
				parameters->system[0].bh[1].theta);
		fprintf(file, text, parameters->system[0].bh[1].varphi, parameters->system[0].coaPhase,
				parameters->system[0].coaTime);
		fprintf(file, "%s %d %s\n", parameters->phase[0], parameters->amp_Code, parameters->spin[0]);
		fprintf(file, text, parameters->system[1].M, parameters->system[1].eta,
				parameters->system[1].chirpM);
		fprintf(file, text, parameters->system[1].mu, parameters->system[1].bh[0].m,
				parameters->system[1].bh[1].m);
		fprintf(file, text, parameters->system[1].incl, parameters->system[1].dist,
				parameters->system[1].F.pol);
		fprintf(file, text, parameters->system[1].F.alpha, parameters->system[1].F.dec,
				parameters->system[1].F.gmst);
		fprintf(file, text, parameters->system[1].bh[0].chi_Amp, parameters->system[1].bh[0].kappa,
				parameters->system[1].bh[0].psi);
		fprintf(file, text, parameters->system[1].bh[1].chi_Amp, parameters->system[1].bh[1].kappa,
				parameters->system[1].bh[1].psi);
		fprintf(file, text, parameters->system[1].bh[0].chi[0], parameters->system[1].bh[0].chi[1],
				parameters->system[1].bh[0].chi[2]);
		fprintf(file, text, parameters->system[1].bh[1].chi[0], parameters->system[1].bh[1].chi[1],
				parameters->system[1].bh[1].chi[2]);
		fprintf(file, text, parameters->system[1].bh[0].theta, parameters->system[1].bh[0].varphi,
				parameters->system[1].bh[1].theta);
		fprintf(file, text, parameters->system[1].bh[1].varphi, parameters->system[1].coaPhase,
				parameters->system[1].coaTime);
		fprintf(file, "%s %d %s ", parameters->phase[1], parameters->amp_Code, parameters->spin[1]);
		fprintf(file, text, parameters->match_Typ, parameters->match_Best,
				parameters->match_Minimax);
		fprintf(file, "%d\n", index);
		fclose(file);
	} else {
		printf("Can not open file: %s, terminating!!!\n", file_Name);
		fflush(stdout);
		abort();
	}
	if (is_First) {
		sprintf(file_Name, "%s/wave%d.txtf", prog->folder, index);
	} else {
		sprintf(file_Name, "%s/wave%d.txt", prog->folder, index);
	}
	file = fopen(file_Name, "w");
	sprintf(temp, "%%%d.%dlg\t", prog->width_Of_Number_To_Plot, prog->precision_To_Plot);
	sprintf(text, "%s%s%s", temp, temp, temp);
	fprintf(file, "#");
	fprintf(file, text, parameters->system[0].M, parameters->system[0].bh[0].m
			/ parameters->system[0].bh[1].m, parameters->system[0].eta);
	fprintf(file, text, parameters->system[0].bh[0].chi_Amp, parameters->system[0].bh[0].kappa,
			parameters->system[0].bh[0].psi);
	fprintf(file, text, parameters->system[0].bh[1].chi_Amp, parameters->system[0].bh[1].kappa,
			parameters->system[0].bh[1].psi);
	fprintf(file, "%s %d %s\n", parameters->phase[0], parameters->amp_Code, parameters->spin[0]);
	fprintf(file, "#");
	fprintf(file, text, parameters->system[1].M, parameters->system[1].bh[0].m
			/ parameters->system[1].bh[1].m, parameters->system[1].eta);
	fprintf(file, text, parameters->system[1].bh[0].chi_Amp, parameters->system[1].bh[0].kappa,
			parameters->system[1].bh[0].psi);
	fprintf(file, text, parameters->system[1].bh[1].chi_Amp, parameters->system[1].bh[1].kappa,
			parameters->system[1].bh[1].psi);
	fprintf(file, "%s %d %s", parameters->phase[1], parameters->amp_Code, parameters->spin[1]);
	fprintf(file, text, parameters->match_Typ, parameters->match_Best, parameters->match_Minimax);
	fprintf(file, "\n");
	long i;
	for (i = 0; i < parameters->min_Length; i++) {
		fprintf(file, "%*.*lg\t", prog->width_Of_Number_To_Plot, prog->precision_To_Plot,
				(double) i * parameters->time_Sampling);
		fprintf(file, text, sig->signal[H1P][i], sig->signal[H1C][i], sig->signal[H1P][i]
				* parameters->system[0].F.F[0] + sig->signal[H1C][i] * parameters->system[0].F.F[1]);
		fprintf(file, text, sig->signal[H2P][i], sig->signal[H2C][i], sig->signal[H2P][i]
				* parameters->system[0].F.F[0] + sig->signal[H2C][i] * parameters->system[0].F.F[1]);
		fprintf(file, "\n");
	}
	if (parameters->shorter) {
		for (; i < parameters->max_Length; i++) {
			fprintf(file, "%*.*lg\t", prog->width_Of_Number_To_Plot, prog->precision_To_Plot,
					(double) i * parameters->time_Sampling);
			fprintf(file, text, sig->signal[H1P][i], sig->signal[H1C][i], sig->signal[H1P][i]
					* parameters->system[0].F.F[0] + sig->signal[H1C][i]
					* parameters->system[0].F.F[1]);
			fprintf(file, "\n");
		}
	} else {
		fprintf(file, "%*.*lg\t", prog->width_Of_Number_To_Plot, prog->precision_To_Plot,
				(double) i * parameters->time_Sampling);
		fprintf(file, "%*s\t%*s\t%*s\t", prog->width_Of_Number_To_Plot, "",
				prog->width_Of_Number_To_Plot, "", prog->width_Of_Number_To_Plot, "");
		fprintf(file, text, sig->signal[H2P][i], sig->signal[H2C][i], sig->signal[H2P][i]
				* parameters->system[0].F.F[0] + sig->signal[H2C][i] * parameters->system[0].F.F[1]);
		fprintf(file, "\n");
	}
	fclose(file);
}

void write_Params_To_File(Program_Parameters *prog, System_Parameters *parameters, char *file_Name) {
	assert(prog);
	assert(parameters);
	static char temp[FILE_NAME_LENGTH];
	static char text[FILE_NAME_LENGTH];
	sprintf(temp, "%%%d.%dlg\t", prog->width_Of_Number, prog->precision);
	sprintf(text, "%s%s%s", temp, temp, temp);
	FILE *file = sfopen(file_Name, "w");
	if (file) {
		fprintf(file, text, parameters->system[0].M, parameters->system[0].eta,
				parameters->system[0].chirpM);
		fprintf(file, text, parameters->system[0].mu, parameters->system[0].bh[0].m,
				parameters->system[0].bh[1].m);
		fprintf(file, text, parameters->system[0].incl, parameters->system[0].dist,
				parameters->system[0].F.pol);
		fprintf(file, text, parameters->system[0].F.alpha, parameters->system[0].F.dec,
				parameters->system[0].F.gmst);
		fprintf(file, text, parameters->system[0].bh[0].chi_Amp, parameters->system[0].bh[0].kappa,
				parameters->system[0].bh[0].psi);
		fprintf(file, text, parameters->system[0].bh[1].chi_Amp, parameters->system[0].bh[1].kappa,
				parameters->system[0].bh[1].psi);
		fprintf(file, text, parameters->system[0].bh[0].chi[0], parameters->system[0].bh[0].chi[1],
				parameters->system[0].bh[0].chi[2]);
		fprintf(file, text, parameters->system[0].bh[1].chi[0], parameters->system[0].bh[1].chi[1],
				parameters->system[0].bh[1].chi[2]);
		fprintf(file, text, parameters->system[0].bh[0].theta, parameters->system[0].bh[0].varphi,
				parameters->system[0].bh[1].theta);
		fprintf(file, text, parameters->system[0].bh[1].varphi, parameters->system[0].coaPhase,
				parameters->system[0].coaTime);
		fprintf(file, "%s %d %s\n", parameters->phase[0], parameters->amp_Code, parameters->spin[0]);
		fprintf(file, text, parameters->system[1].M, parameters->system[1].eta,
				parameters->system[1].chirpM);
		fprintf(file, text, parameters->system[1].mu, parameters->system[1].bh[0].m,
				parameters->system[1].bh[1].m);
		fprintf(file, text, parameters->system[1].incl, parameters->system[1].dist,
				parameters->system[1].F.pol);
		fprintf(file, text, parameters->system[1].F.alpha, parameters->system[1].F.dec,
				parameters->system[1].F.gmst);
		fprintf(file, text, parameters->system[1].bh[0].chi_Amp, parameters->system[1].bh[0].kappa,
				parameters->system[1].bh[0].psi);
		fprintf(file, text, parameters->system[1].bh[1].chi_Amp, parameters->system[1].bh[1].kappa,
				parameters->system[1].bh[1].psi);
		fprintf(file, text, parameters->system[1].bh[0].chi[0], parameters->system[1].bh[0].chi[1],
				parameters->system[1].bh[0].chi[2]);
		fprintf(file, text, parameters->system[1].bh[1].chi[0], parameters->system[1].bh[1].chi[1],
				parameters->system[1].bh[1].chi[2]);
		fprintf(file, text, parameters->system[1].bh[0].theta, parameters->system[1].bh[0].varphi,
				parameters->system[1].bh[1].theta);
		fprintf(file, text, parameters->system[1].bh[1].varphi, parameters->system[1].coaPhase,
				parameters->system[1].coaTime);
		fprintf(file, "%s %hd %s ", parameters->phase[1], parameters->amp_Code, parameters->spin[1]);
		fprintf(file, text, parameters->match_Typ, parameters->match_Best,
				parameters->match_Minimax);
		fclose(file);
	} else {
		printf("Can not open file: %s, terminating!!!\n", file_Name);
		fflush(stdout);
		abort();
	}
}
