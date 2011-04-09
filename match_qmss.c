/**
 * @file match_qmss.c
 *
 * @date Feb 25, 2011
 * @author vereb
 */

#include "match_qmss.h"

void run_Algorithm(Program_Parameters *program_Parameters, System_Parameters *parameters,
		binary_System *limits) {
	assert(program_Parameters);
	assert(parameters);
	assert(limits);
	assert(program_Parameters->number_Of_Runs > 0);
	srand(86);
	short is_Good;
	for (short i = 0; i < program_Parameters->number_Of_Runs;) {
		generate_Parameters(parameters, limits);
		is_Good = incrementing_Spins(program_Parameters, parameters, i);
		if (is_Good) {
			i++;
		}
		if ((i + 1) % 10 == 0) {
			printf("%ld%%\n", 100 * i / program_Parameters->number_Of_Runs);
		}
	}
}

void generate_Same_Parameters(System_Parameters *parameters, binary_System *limits) {
	assert(parameters);
	assert(limits);
	gen_Parameters(&parameters->system[0], &limits[0], &limits[1], ETAM, KAPPA_PSI);
	memcpy(&parameters->system[1], &parameters->system[0], sizeof(binary_System));
}

void generate_Parameters(System_Parameters *parameters, binary_System *limits) {
	assert(parameters);
	assert(limits);
	gen_Parameters(&parameters->system[0], &limits[0], &limits[1], ETAM, KAPPA_PSI);
	gen_Parameters(&parameters->system[1], &limits[0], &limits[1], ETAM, KAPPA_PSI);
}

short incrementing_Spins(Program_Parameters *prog, System_Parameters* parameters, short index) {
	assert(prog);
	assert(parameters);
	FILE*file;
	char file_Name[FILENAME_MAX];
	short is_Good;
	parameters->critical_Match = 0.0;
	signalStruct first, second;
	is_Good = calc_Matches_For_ParameterPair(prog, parameters, &first);
	if (!is_Good) {
		return NOT_FOUND;
	}
	parameters->critical_Match = parameters->match_Minimax;
	double min_Spin = parameters->system[MOD_SPIN_INDEX].bh[1].chi_Amp;
	for (; parameters->system[MOD_SPIN_INDEX].bh[0].chi_Amp < parameters->max_Spin; parameters->system[MOD_SPIN_INDEX].bh[0].chi_Amp
			+= parameters->spin_Step) {
		for (; parameters->system[MOD_SPIN_INDEX].bh[1].chi_Amp < parameters->max_Spin; parameters->system[MOD_SPIN_INDEX].bh[1].chi_Amp
				+= parameters->spin_Step) {
			if (parameters->system[MOD_SPIN_INDEX].bh[0].chi_Amp == min_Spin
					&& parameters->system[MOD_SPIN_INDEX].bh[1].chi_Amp == min_Spin) {
				continue;
			}
			is_Good = calc_Matches_For_ParameterPair(prog, parameters, &second);
			if (is_Good && parameters->match_Minimax > parameters->critical_Match) {
				sprintf(file_Name, "%s/first%hd.txt", prog->folder, index);
				file = safely_Open_File_For_Writing(file_Name);
				print_Two_Signals(file, &first, parameters->time_Sampling,
						prog->width_Of_Number_To_Plot, prog->precision_To_Plot);
				fclose(file);
				sprintf(file_Name, "%s/found%hd.txt", prog->folder, index);
				file = safely_Open_File_For_Writing(file_Name);
				print_Two_Signals(file, &first, parameters->time_Sampling,
						prog->width_Of_Number_To_Plot, prog->precision_To_Plot);
				fclose(file);
				destroy_Signal_Struct(&first);
				destroy_Signal_Struct(&second);
				return FOUND;
			}
			destroy_Signal_Struct(&second);
		}
		parameters->system[MOD_SPIN_INDEX].bh[1].chi_Amp = min_Spin;
	}
	destroy_Signal_Struct(&first);
	return NOT_FOUND;
}

inline void increment_Spin_Of_Binary_System(binary_System *system, double step) {
	assert(system);
	assert(step > 0.);
	system->bh[0].chi_Amp += step;
	system->bh[1].chi_Amp += step;
}

inline void increment_Spins(System_Parameters* parameters) {
	assert(parameters);
	increment_Spin_Of_Binary_System(&parameters->system[0], parameters->spin_Step);
	increment_Spin_Of_Binary_System(&parameters->system[1], parameters->spin_Step);
}

short calc_Matches_For_ParameterPair(Program_Parameters *prog, System_Parameters *parameters,
		signalStruct *sig) {
	assert(prog);
	assert(parameters);
	static LALParameters lalparams;
	initLALParameters(&lalparams, parameters);
	for (short i = 0; i < 2; i++) {
		memset(&lalparams.status, 0, sizeof(LALStatus));
		LALGenerateInspiral(&lalparams.status, &lalparams.waveform[i], &lalparams.injParams[i],
				&lalparams.ppnParams);
		if (lalparams.status.statusCode) {
			fprintf(stderr, "%d: LALSQTPNWaveformTest: error generating waveform\n", i);
			XLALSQTPNDestroyCoherentGW(&lalparams.waveform[0]);
			return NOT_FOUND;
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
	parameters->freq_Step = 1. / (lalparams.ppnParams.deltaT * lalparams.max_Length);
	create_Signal_Struct(sig, lalparams.waveform[!lalparams.shorter].f->data->length);
	createPSD(&lalparams, sig);
	for (short i = 0; i < 2; i++) {
		for (long j = 0; j < lalparams.waveform[i].f->data->length; j++) {
			set_Signal_From_A1A2(i, j, sig, &lalparams,
					parameters->system[i].F.antenna_Beam_Pattern);
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
	if (minfr == maxfr) {
		fprintf(stderr, "The two frequencies are too close!!! %lg~%lg: %ld %ld\n",
				parameters->freq_Min, parameters->freq_Max, minfr, maxfr);
		XLALSQTPNDestroyCoherentGW(&lalparams.waveform[0]);
		XLALSQTPNDestroyCoherentGW(&lalparams.waveform[1]);
		return NOT_FOUND;
	}
	parameters->match_Typ = parameters->match_Best = parameters->match_Minimax = 0.0;
	calc_Matches(sig, minfr, maxfr, &parameters->match_Typ, &parameters->match_Best,
			&parameters->match_Minimax);
	XLALSQTPNDestroyCoherentGW(&lalparams.waveform[0]);
	XLALSQTPNDestroyCoherentGW(&lalparams.waveform[1]);
	XLALFree(lalparams.randIn.psd.data);
	if (parameters->match_Minimax < parameters->min_Match) {
		return NOT_FOUND;
	}
	return FOUND;
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
		lalparams->injParams[i].amp_order = parameters->amp_Code[i];
		snprintf(lalparams->injParams[i].waveform, LIGOMETA_WAVEFORM_MAX * sizeof(CHAR), "%s%s%s",
				parameters->approx[i], parameters->phase[i], parameters->spin[i]);
	}
}

void createPSD(LALParameters *lalparams, signalStruct *sig) {
	assert(lalparams);
	assert(sig);
	lalparams->randIn.psd.length = lalparams->max_Length;
	double df = 1. / lalparams->ppnParams.deltaT / lalparams->randIn.psd.length;
	lalparams->randIn.psd.data = (REAL8*)LALMalloc(sizeof(REAL8) * lalparams->randIn.psd.length);
	LALNoiseSpectralDensity(&lalparams->status, &lalparams->randIn.psd, &LALLIGOIPsd, df);
	for (long j = 0; j < lalparams->randIn.psd.length; j++) {
		sig->psd[j] = lalparams->randIn.psd.data[j];
	}
}

void set_Signal_From_A1A2(short index, long elem, signalStruct *sig, LALParameters *lal, double F[]) {
	double a1, a2, phi, shift;
	a1 = lal->waveform[index].a->data->data[2 * elem];
	a2 = lal->waveform[index].a->data->data[2 * elem + 1];
	phi = lal->waveform[index].phi->data->data[elem] - lal->waveform[index].phi->data->data[0];
	shift = 0.0;
	sig->signal[2 * index][elem] = a1 * M_SQRT2 / 2.0;
	sig->signal[2 * index + 1][elem] = a2 * M_SQRT2 / 2.0;
	sig->signal[RESPONSE1 + index][elem] = sig->signal[2 * index][elem] * F[0] + sig->signal[2
			* index + 1][elem] * F[1];
}

void set_Signal_FromHPHC(short index, long elem, signalStruct *sig, LALParameters *lal, double F[]) {
	sig->signal[2 * index][elem] = lal->waveform[index].h->data->data[2 * elem];
	sig->signal[2 * index + 1][elem] = lal->waveform[index].h->data->data[2 * elem + 1];
	sig->signal[RESPONSE1 + index][elem] = sig->signal[2 * index][elem] * F[0] + sig->signal[2
			* index + 1][elem] * F[1];
}

void readExactParameters(FILE *file, System_Parameters *params) {
	char name[100];
	for (short i = 0; i < 2; i++) {
		fscanf(file, "%s ", name);
		fscanf(file, "%lg ", &params->system[i].bh[0].m);
		fscanf(file, "%lg ", &params->system[i].bh[1].m);
		fscanf(file, "%lg ", &params->system[i].bh[0].chi_Amp);
		fscanf(file, "%lg ", &params->system[i].bh[0].kappa);
		params->system[i].bh[0].kappa *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
		fscanf(file, "%lg ", &params->system[i].bh[0].psi);
		params->system[i].bh[0].psi *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
		fscanf(file, "%lg ", &params->system[i].bh[1].chi_Amp);
		fscanf(file, "%lg ", &params->system[i].bh[1].kappa);
		params->system[i].bh[1].kappa *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
		fscanf(file, "%lg ", &params->system[i].bh[1].psi);
		params->system[i].bh[1].psi *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
		fscanf(file, "%lg ", &params->system[i].dist);
		fscanf(file, "%lg ", &params->system[i].incl);
		params->system[i].incl *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
		fscanf(file, "%lg ", &params->freq_Sampling);
		params->time_Sampling = 1.0 / params->freq_Sampling;
		fscanf(file, "%lg ", &params->freq_Initial);
		fscanf(file, "%lg ", &params->freq_Max);
		fscanf(file, "%s ", params->phase[i]);
		fscanf(file, "%s ", params->spin[i]);
		fscanf(file, "%hd ", &params->amp_Code[i]);
		fscanf(file, "%s\n", params->approx[i]);
		params->system[i].F.declination = params->system[i].F.gmst
				= params->system[i].F.greenwich_Hour_Angle = params->system[i].F.polarization
						= params->system[i].F.right_Ascention = 0.;
		convert_Masses(&params->system[i], FROM_M1M2);
		convert_Spins(&params->system[i], FROM_KAPPA_PSI);
		calc_Antenna_Pattern_For(LH, &params->system[i].F);
	}
}

void read_Program_Parameters(Program_Parameters *parameters, System_Parameters *params,
		char *file_Name) {
	assert(parameters);
	assert(params);
	assert(file_Name);
	FILE *file = fopen(file_Name, "r");
	fscanf(file, "%ld\n", &parameters->number_Of_Runs);
	fscanf(file, "%hd\n", &parameters->precision);
	parameters->width_Of_Number = parameters->precision + SPECIAL_CHARACTER_LENGTH;
	fscanf(file, "%hd\n", &parameters->precision_To_Plot);
	parameters->width_Of_Number_To_Plot = parameters->precision_To_Plot + SPECIAL_CHARACTER_LENGTH;
	fscanf(file, "%s\n", parameters->folder);
	fscanf(file, "%lg\n", &params->min_Match);
	fscanf(file, "%lg\n", &params->max_Spin);
	fscanf(file, "%lg\n", &params->spin_Step);
	fscanf(file, "%lg\n", &params->freq_Sampling);
	params->time_Sampling = 1. / params->freq_Sampling;
	fscanf(file, "%lg\n", &params->freq_Initial);
	fscanf(file, "%lg\n", &params->freq_Max);
	fscanf(file, "%lg\n", &params->delta_Length);
	fscanf(file, "%s\n", params->approx[0]);
	fscanf(file, "%s\n", params->phase[0]);
	fscanf(file, "%hd\n", &params->amp_Code[0]);
	fscanf(file, "%s\n", params->spin[0]);
	fscanf(file, "%s\n", params->approx[1]);
	fscanf(file, "%s\n", params->phase[1]);
	fscanf(file, "%hd\n", &params->amp_Code[1]);
	fscanf(file, "%s\n", params->spin[1]);
	fclose(file);
}
