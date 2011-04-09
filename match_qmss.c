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
	for (; parameters->system[MOD_SPIN_INDEX].bh[0].chi_Amp < prog->max_Spin; parameters->system[MOD_SPIN_INDEX].bh[0].chi_Amp
			+= prog->spin_Step) {
		for (; parameters->system[MOD_SPIN_INDEX].bh[1].chi_Amp < prog->max_Spin; parameters->system[MOD_SPIN_INDEX].bh[1].chi_Amp
				+= prog->spin_Step) {
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

inline void increment_Spins(System_Parameters* parameters, double step) {
	assert(parameters);
	increment_Spin_Of_Binary_System(&parameters->system[0], step);
	increment_Spin_Of_Binary_System(&parameters->system[1], step);
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
	createPSD(&lalparams, sig->psd);
	for (short i = 0; i < 2; i++) {
		setSignal_From_A1A2(i, sig, &lalparams, parameters->system[i].F.antenna_Beam_Pattern);
	}
	double fr = 0.;
	long minfr = 0, maxfr = 0;
	while (fr < parameters->freq_Min) {
		fr += parameters->freq_Step;
		maxfr = ++minfr;
	}
	while (fr < prog->freq_Max) {
		fr += parameters->freq_Step;
		maxfr++;
	}
	if (minfr == maxfr) {
		fprintf(stderr, "The two frequencies are too close!!! %lg~%lg: %ld %ld\n",
				parameters->freq_Min, prog->freq_Max, minfr, maxfr);
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
	if (parameters->match_Minimax < prog->min_Match) {
		return NOT_FOUND;
	}
	return FOUND;
}
