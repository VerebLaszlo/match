/**
 * @file parameters.c
 *
 * @date Apr 9, 2011
 * @author vereb
 */

#include <assert.h>
#include "parameters.h"

void readExactParameters(FILE *file, SystemParameter *params) {
	char stringFormat[] = "%25s ";
	for (ushort i = 0; i < NUMBER_OF_SYSTEMS; i++) {
		fscanf(file, stringFormat, params->name[i]);
		fscanf(file, "%lg ", &params->system[i].mass.mass[0]);
		fscanf(file, "%lg ", &params->system[i].mass.mass[1]);
		fscanf(file, "%lg ", &params->system[i].spin[0].magnitude);
		fscanf(file, "%lg ", &params->system[i].spin[0].inclination);
		params->system[i].spin[0].inclination = radianFromDegree(
			params->system[i].spin[0].inclination);
		fscanf(file, "%lg ", &params->system[i].spin[0].azimuth);
		params->system[i].spin[0].azimuth = radianFromDegree(params->system[i].spin[0].azimuth);
		fscanf(file, "%lg ", &params->system[i].spin[1].magnitude);
		fscanf(file, "%lg ", &params->system[i].spin[1].inclination);
		params->system[i].spin[1].inclination = radianFromDegree(
			params->system[i].spin[1].inclination);
		fscanf(file, "%lg ", &params->system[i].spin[1].azimuth);
		params->system[i].spin[1].azimuth = radianFromDegree(params->system[i].spin[1].azimuth);
		fscanf(file, "%lg ", &params->system[i].distance);
		fscanf(file, "%lg ", &params->system[i].inclination);
		params->system[i].inclination = radianFromDegree(params->system[i].inclination);
		fscanf(file, "%lg ", &params->samplingFrequency);
		params->samplingTime = 1.0 / params->samplingFrequency;
		fscanf(file, "%lg ", &params->initialFrequency);
		fscanf(file, stringFormat, params->approximant[i]);
		fscanf(file, stringFormat, params->phase[i]);
		fscanf(file, stringFormat, params->spin[i]);
		fscanf(file, stringFormat, &params->amplitude[i]);
		params->detector[i].declination = params->detector[i].greenwichMeanSiderealTime = params
			->detector[i].greenwichHourAngle = params->detector[i].polarization =
			params->detector[i].rightAscention = 0.0;
	}
}

void readSystemParameters(FILE *file, SystemParameter *params) {
	short length = 100;
	char line[length];
	fgets(line, length, file);
	sscanf(line, "%lg\n", &params->initialFrequency);
	fgets(line, length, file);
	sscanf(line, "%lg\n", &params->samplingFrequency);
	params->samplingTime = 1. / params->samplingFrequency;
	fgets(line, length, file);
	fgets(line, length, file);
	fgets(line, length, file);
	sscanf(line, "%25s %25s\n", params->approximant[0], params->approximant[1]);
	fgets(line, length, file);
	sscanf(line, "%25s %25s %*s\n", params->phase[0], params->phase[1]);
	fgets(line, length, file);
	sscanf(line, "%25s %25s %*s\n", params->spin[0], params->spin[1]);
	fgets(line, length, file);
	sscanf(line, "%25s %25s %*s\n", &params->amplitude[0], &params->amplitude[1]);
	fgets(line, length, file);
	read_Binary_Parameter_Limits(file, params->system);
}

void readProgramParameters(FILE *file, ProgramParameter *params) {
	assert(file);
	assert(params);
	short length = 100;
	char line[length];
	fgets(line, length, file);
	sscanf(line, "%ld\n", &params->numberOfRuns);
	fgets(line, length, file);
	sscanf(line, "%hd\n", &params->precision[TO_PLOT]);
	fgets(line, length, file);
	sscanf(line, "%hd\n", &params->width[TO_PLOT]);
	fgets(line, length, file);
	sscanf(line, "%hd\n", &params->precision[TO_BACKUP]);
	fgets(line, length, file);
	sscanf(line, "%hd\n", &params->width[TO_BACKUP]);
	fgets(line, length, file);
	sscanf(line, "%4096s\n", params->outputDirectory);
}

void print_Program_Parameters(FILE*file, ProgramParameters *params) {
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

void print_System_Parameters(FILE *file, System_Parameters *params) {
	fprintf(file, "%10s %10.4lg\n", "freq_I", params->freq_Initial);
	fprintf(file, "%10s %10.4lg\n", "freq_S", params->freq_Sampling);
	fprintf(file, "%10s %10s %10s\n", "approx", params->approx[0], params->approx[1]);
	fprintf(file, "%10s %10s %10s\n", "phase", params->phase[0], params->phase[1]);
	fprintf(file, "%10s %10s %10s\n", "spin", params->spin[0], params->spin[1]);
	fprintf(file, "%10s %10d %10d\n", "amp", params->amp_Code[0], params->amp_Code[1]);
	print_Binary_Parameter_Limits(file, params->system);
}

void print_System_Parameters_For_Plot(FILE *file, System_Parameters *params) {
	print_Binary_Parameter_For_Plot(file, params->system);
	fprintf(file, "%-13s %10.4lg %10.4lg %10.4lg\n", "#matches    ", params->match_Typ,
		params->match_Minimax, params->match_Best);
}
