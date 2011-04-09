/**
 * @file parameters.c
 *
 * @date Apr 9, 2011
 * @author vereb
 */

#include "parameters.h"

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

void read_Program_Parameters(FILE*file, Program_Parameters *params) {
	assert(file);
	assert(params);
	short length = 100;
	char line[length];
	fgets(line, length, file);
	sscanf(line, "%ld\n", &params->number_Of_Runs);
	fgets(line, length, file);
	sscanf(line, "%hd\n", &params->precision);
	params->width_Of_Number = params->precision + SPECIAL_CHARACTER_LENGTH;
	fgets(line, length, file);
	sscanf(line, "%hd\n", &params->precision_To_Plot);
	params->width_Of_Number_To_Plot = params->precision_To_Plot + SPECIAL_CHARACTER_LENGTH;
	fgets(line, length, file);
	sscanf(line, "%s\n", params->folder);
	fgets(line, length, file);
	sscanf(line, "%lg\n", &params->min_Match);
	fgets(line, length, file);
	sscanf(line, "%lg\n", &params->max_Spin);
	fgets(line, length, file);
	sscanf(line, "%lg\n", &params->spin_Step);
	fgets(line, length, file);
	sscanf(line, "%lg\n", &params->freq_Max);
	fgets(line, length, file);
	sscanf(line, "%lg\n", &params->delta_Length);
}

void print_Program_Parameters(FILE*file, Program_Parameters *params) {
	puts("T");
	fprintf(file, "%10s %10ld\n", "numOfRuns", params->number_Of_Runs);
	fprintf(file, "%10s %10hd\n", "prec", params->precision);
	fprintf(file, "%10s %10d\n", "width", params->width_Of_Number);
	fprintf(file, "%10s %10hd\n", "precPlot", params->precision_To_Plot);
	fprintf(file, "%10s %10d\n", "widthPlot", params->width_Of_Number_To_Plot);
	fprintf(file, "%10s %10s\n", "folder", params->folder);
	fprintf(file, "%10s %10.4lg\n", "min match", params->min_Match);
	fprintf(file, "%10s %10.4lg\n", "max spin", params->max_Spin);
	fprintf(file, "%10s %10.4lg\n", "spin step", params->spin_Step);
	fprintf(file, "%10s %10.4lg\n", "freq max", params->freq_Max);
	fprintf(file, "%10s %10.4lg\n", "delta L", params->delta_Length);
}

void read_System_Parameters(FILE *file, System_Parameters *params) {
	short length = 100;
	char line[length];
	fgets(line, length, file);
	sscanf(line, "%lg\n", &params->freq_Initial);
	fgets(line, length, file);
	sscanf(line, "%lg\n", &params->freq_Sampling);
	params->time_Sampling = 1. / params->freq_Sampling;
	fgets(line, length, file);
	fgets(line, length, file);
	fgets(line, length, file);
	sscanf(line, "%s %s %*s\n", params->approx[0], params->approx[1]);
	fgets(line, length, file);
	sscanf(line, "%s %s %*s\n", params->phase[0], params->phase[1]);
	fgets(line, length, file);
	sscanf(line, "%s %s %*s\n", params->spin[0], params->spin[1]);
	fgets(line, length, file);
	sscanf(line, "%hd %hd %*s\n", &params->amp_Code[0], &params->amp_Code[1]);
	fgets(line, length, file);
	read_Binary_Parameter_Limits(file, params->system);
}

void print_System_Parameters(FILE *file, System_Parameters *params) {
	fprintf(file, "%10s %10.4lg\n", "freq_I", params->freq_Initial);
	fprintf(file, "%10s %10.4lg\n", "freq_S", params->freq_Sampling);
	fprintf(file, "%10s %10s %10s\n", "approx", params->approx[0], params->approx[1]);
	fprintf(file, "%10s %10s %10s\n", "phase", params->phase[0], params->phase[1]);
	fprintf(file, "%10s %10s %10s\n", "spin", params->spin[0], params->spin[1]);
	fprintf(file, "%10s %10d %10d\n", "amp", params->amp_Code[0], params->amp_Code[1]);
	print_Binary_Parameter_Limits(file, params->system);
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
