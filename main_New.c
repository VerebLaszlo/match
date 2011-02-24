/**
 * @file main_New.c
 *
 * @date Feb 17, 2011
 * @author vereb
 */

#include "match.h"

#define MOD_SPIN_INDEX 0///<C
#include <lal/LALSQTPNWaveformInterface.h>
#include <lal/LALNoiseModelsInspiral.h>
//#include <lal/GeneratePPNInspiral.h>
#include <lal/SimulateCoherentGW.h>	// CoherentGW
//#include <lal/LIGOMetadataTables.h>		// SimInspiralTable
#include <lal/GenerateInspiral.h>
#include <lal/LALInspiralBank.h>
#include <lal/LALDatatypes.h>		// LALStatus)
#include <lal/LALInspiral.h>
#include <lal/LALStdlib.h>
#include <lal/RealFFT.h>

short is_First;///<a

/**
 * X
 */
typedef enum Constants {
	EXTRA_CHARACTERS = 6, FILE_NAME_LENGTH = 100,
} Constants;

/**
 * X
 */
typedef struct Program_Parameters {
	char (*output_Directories)[FILE_NAME_LENGTH];///<a
	long number_Of_Runs;///<a
	short precision;///<a
	short precision_To_Plot;///<a
	short width_Of_Number;///<a
	short width_Of_Number_To_Plot;///<a
	char folder[FILE_NAME_LENGTH];///<a
} Program_Parameters;

/**
 * X
 */
typedef struct System_Parameters {
	binary_System system[2];///<a
	double max_Spin;///<a
	double spin_Step;///<a
	double freq_Sampling;///<a
	double freq_Initial;///<a
	double time_Sampling;///<a
	double match_Typ;///<aa
	double match_Best;///<a
	double match_Minimax;///<a
	short shorter;///<a
	long min_Length;///<a
	long max_Length;///<a
	double freq_Min;///<a
	double freq_Max;///<a
	double freq_Step;///<a
	double min_Match;///<a
	double critical_Match;///<a
	double delta_Length;///<a
	char approx[2][FILE_NAME_LENGTH];///<a
	char phase[2][FILE_NAME_LENGTH];///<a
	char amp[2][FILE_NAME_LENGTH];///<a
	char spin[2][FILE_NAME_LENGTH];///<a
} System_Parameters;

/**
 * X
 */
typedef struct LALParameters {
	LALStatus status;///<a
	CoherentGW waveform[2];///<a
	SimInspiralTable injParams[2];///<a
	PPNParamStruc ppnParams;///<a
	RandomInspiralSignalIn randIn;///<a
	short shorter;///<a
	long min_Length;///<a
	long max_Length;///<a
} LALParameters;

void read_Program_Parameters(Program_Parameters *parameters, System_Parameters *params,
		char *file_Name);

void read_Parameters(binary_System *parameters, char *file_Name);

void run_Algorithm(Program_Parameters *program_Parameters, System_Parameters *parameters,
		binary_System *limits);

void generate_Parameters(System_Parameters *parameters, binary_System *limits);

short incrementing_Spins(Program_Parameters *prog, System_Parameters* parameters);

void increment_Spin_Of_Binary_System(binary_System *system, double step);

void increment_Spins(System_Parameters* parameters);

short calc_Matches_For_ParameterPair(Program_Parameters *prog, System_Parameters *parameters);

void initLALParameters(LALParameters *lalparams, System_Parameters *parameters);

void createPSD(LALParameters *lalparams, signalStruct *sig);

void write_Waves_To_Files(Program_Parameters *prog, System_Parameters *parameters,
		signalStruct *sig);

void write_Wave_To_File(Program_Parameters *prog, System_Parameters *parameters, signalStruct *sig,
		short index);

/**
 * Done.
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char *argv[]) {
	char program_Parameters_File_Name[FILE_NAME_LENGTH];
	char parameters_File_Name[FILE_NAME_LENGTH];
	if (argc != 3) {
		puts("\"file for program parameters\" \"file for the parameter limits\"!!!");
		exit(-1);
	}
	sprintf(program_Parameters_File_Name, argv[1]);
	sprintf(parameters_File_Name, argv[2]);
	Program_Parameters program_Parameters;
	binary_System limits_Of_Parameters[2];
	System_Parameters parameters;
	read_Program_Parameters(&program_Parameters, &parameters, program_Parameters_File_Name);
	read_Parameters(limits_Of_Parameters, parameters_File_Name);
	puts("Start!!");
	run_Algorithm(&program_Parameters, &parameters, limits_Of_Parameters);
	puts("Done!!!");
	return 0;
}

/**
 * Done.
 * @param parameters
 * @param params
 * @param file_Name
 */
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
	fscanf(file, "%s\n", params->amp[0]);
	fscanf(file, "%s\n", params->spin[0]);
	fscanf(file, "%s\n", params->approx[1]);
	fscanf(file, "%s\n", params->phase[1]);
	fscanf(file, "%s\n", params->amp[1]);
	fscanf(file, "%s\n", params->spin[1]);
	fclose(file);
}

/**
 * Done.
 * @param parameters
 * @param file_Name
 */
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
	fscanf(file, "%lg %lg\n", &parameters[0].bh[1].kappa, &parameters[1].bh[1].kappa);
	fscanf(file, "%lg %lg\n", &parameters[0].bh[0].psi, &parameters[1].bh[0].psi);
	fscanf(file, "%lg %lg\n", &parameters[0].bh[1].psi, &parameters[1].bh[1].psi);
	fscanf(file, "%lg %lg\n", &parameters[0].incl, &parameters[1].incl);
	fscanf(file, "%lg %lg\n", &parameters[0].dist, &parameters[1].dist);
	fscanf(file, "%lg %lg\n", &parameters[0].F.pol, &parameters[1].F.pol);
	fscanf(file, "%lg %lg\n", &parameters[0].F.alpha, &parameters[1].F.alpha);
	fscanf(file, "%lg %lg\n", &parameters[0].F.dec, &parameters[1].F.dec);
	fscanf(file, "%lg %lg\n", &parameters[0].F.gmst, &parameters[1].F.gmst);
	fclose(file);
}

/**
 * Done.
 * @param program_Parameters
 * @param parameters
 * @param limits
 */
void run_Algorithm(Program_Parameters *program_Parameters, System_Parameters *parameters,
		binary_System *limits) {
	assert(program_Parameters);
	assert(parameters);
	assert(limits);
	assert(program_Parameters->number_Of_Runs > 0);
	srand(86);
	short is_Good;
	for (long i = 0; i < program_Parameters->number_Of_Runs;) {
		generate_Parameters(parameters, limits);
		is_Good = incrementing_Spins(program_Parameters, parameters);
		if (is_Good) {
			i++;
		}
	}
}

/**
 * Done
 * @param parameters
 * @param limits
 */
void generate_Parameters(System_Parameters *parameters, binary_System *limits) {
	assert(parameters);
	assert(limits);
	gen_Parameters(&parameters->system[0], &limits[0], &limits[1], ETAM, KAPPA_PSI);
	memcpy(&parameters->system[1], &parameters->system[0], sizeof(binary_System));
}

/**
 * Done.
 * @param prog
 * @param parameters
 * @return
 */
short incrementing_Spins(Program_Parameters *prog, System_Parameters* parameters) {
	assert(prog);
	assert(parameters);
	char temp[FILE_NAME_LENGTH];
	sprintf(temp, "%s", prog->folder);
	is_First = 1;
	short is_Good;
	parameters->critical_Match = 0.0;
	for (; parameters->system[MOD_SPIN_INDEX].bh[0].chi_Amp < parameters->max_Spin; increment_Spin_Of_Binary_System(
			&parameters->system[MOD_SPIN_INDEX], parameters->spin_Step)) {
		is_Good = calc_Matches_For_ParameterPair(prog, parameters);
		break;
		if (is_First) {
			if (!is_Good) {
				return is_Good;
			}
			parameters->critical_Match = parameters->match_Minimax;
			is_First = 0;
		}
		sprintf(prog->folder, "%s", temp);
	}
	return is_Good;
}

/** Done
 * @param system
 * @param step
 */
inline void increment_Spin_Of_Binary_System(binary_System *system, double step) {
	assert(system);
	assert(step>0.);
	system->bh[0].chi_Amp += step;
	system->bh[1].chi_Amp += step;
}

/** Done
 * @param parameters
 */
inline void increment_Spins(System_Parameters* parameters) {
	assert(parameters);
	increment_Spin_Of_Binary_System(&parameters->system[0], parameters->spin_Step);
	increment_Spin_Of_Binary_System(&parameters->system[1], parameters->spin_Step);
}

/**
 * Done.
 * @param prog
 * @param parameters
 * @return
 */
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
			continue;
		}
		parameters->system[i].coaPhase
				= lalparams.waveform[i].phi->data->data[lalparams.waveform[i].phi->data->length];
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
	printf("%ld\n", sig.size);
	printf("%ld, %lg\n", minfr, fr);
	while (fr < parameters->freq_Max) {
		fr += parameters->freq_Step;
		maxfr++;
	}
	printf("%ld, %lg\n", maxfr, fr);
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

/**
 * Done
 * @param lalparams
 * @param parameters
 */
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
		snprintf(lalparams->injParams[i].waveform, LIGOMETA_WAVEFORM_MAX * sizeof(CHAR),
				"%s%s%s%s", parameters->approx[i], parameters->phase[i], parameters->amp[i],
				parameters->spin[i]);
	}
}

/**
 * Done
 * @param lalparams
 * @param sig
 */
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

/**
 * Done
 * @param prog
 * @param parameters
 * @param sig
 */
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
		puts("A");
		sprintf(prog->folder, "%s/best", temp);
		write_Wave_To_File(prog, parameters, sig, index[0]);
		index[0]++;
	} else if (parameters->match_Minimax > parameters->critical_Match) {
		puts("B");
		sprintf(prog->folder, "%s/match", temp);
		write_Wave_To_File(prog, parameters, sig, index[1]);
		index[1]++;
	}
}

/**
 * Done.
 * @todo előbb megnyitni valahol, hogy megkezdjük írással a fájlt, és ne csak hozzáfüzzünk!!!
 * @param prog
 * @param parameters
 * @param sig
 * @param index
 */
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
		fprintf(file, text, parameters->system[0].bh[0].chi_Amp,
				parameters->system[0].bh[1].chi_Amp, parameters->system[1].bh[0].chi_Amp);
		fprintf(file, text, parameters->system[1].bh[1].chi_Amp, parameters->system[0].bh[0].kappa,
				parameters->system[0].bh[0].psi);
		fprintf(file, text, parameters->system[0].bh[0].theta, parameters->system[0].bh[0].varphi,
				parameters->system[0].bh[1].kappa);
		fprintf(file, text, parameters->system[0].bh[1].psi, parameters->system[0].bh[1].theta,
				parameters->system[0].bh[1].varphi);
		fprintf(file, text, parameters->system[0].bh[0].chi[0], parameters->system[0].bh[0].chi[1],
				parameters->system[0].bh[0].chi[2]);
		fprintf(file, text, parameters->system[0].bh[1].chi[0], parameters->system[0].bh[1].chi[1],
				parameters->system[0].bh[1].chi[2]);
		fprintf(file, text, parameters->system[0].coaPhase, parameters->system[1].coaPhase,
				parameters->system[0].coaTime);
		fprintf(file, text, parameters->system[1].coaTime, parameters->system[0].F.F[0],
				parameters->system[0].F.F[0]);
		fprintf(file, "%s %s ", parameters->phase[0], parameters->phase[1]);
		fprintf(file, "%s %s ", parameters->amp[0], parameters->amp[1]);
		fprintf(file, "%8s %8s ", parameters->spin[0], parameters->spin[1]);
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
			/ parameters->system[0].bh[1].m, parameters->system[0].coaTime);
	fprintf(file, text, parameters->system[1].coaTime, parameters->system[1].bh[0].chi_Amp,
			parameters->system[1].bh[1].chi_Amp);
	fprintf(file, text, parameters->system[0].bh[0].chi_Amp, parameters->system[0].bh[1].chi_Amp,
			parameters->system[0].bh[0].kappa);
	fprintf(file, text, parameters->system[0].bh[0].psi, parameters->system[0].bh[1].kappa,
			parameters->system[0].bh[1].psi);
	fprintf(file, "%s %s ", parameters->phase[0], parameters->phase[1]);
	fprintf(file, "%s %s ", parameters->amp[0], parameters->amp[1]);
	fprintf(file, "%8s %8s ", parameters->spin[0], parameters->spin[1]);
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
