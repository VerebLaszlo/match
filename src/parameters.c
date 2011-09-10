/**
 * @file parameters.c
 *
 * @date Apr 9, 2011
 * @author vereb
 */

#include <assert.h>
#include "parameters.h"

void readExactParameters(FILE *file, SystemParameter *params) {
	memset(params, 0, sizeof(SystemParameter));
	char stringFormat[] = "%25s ";
	for (ushort i = 0; i < NUMBER_OF_SYSTEMS; i++) {
		fscanf(file, stringFormat, params->name[i]);
		fscanf(file, "%lg ", &params->system[i].mass.mass[0]);
		fscanf(file, "%lg ", &params->system[i].mass.mass[1]);
		fscanf(file, "%lg ", &params->system[i].spin[0].magnitude);
		fscanf(file, "%lg ", &params->system[i].spin[0].inclination[FIXED]);
		params->system[i].spin[0].inclination[FIXED] = radianFromDegree(
			params->system[i].spin[0].inclination[FIXED]);
		fscanf(file, "%lg ", &params->system[i].spin[0].azimuth[FIXED]);
		params->system[i].spin[0].azimuth[FIXED] = radianFromDegree(
			params->system[i].spin[0].azimuth[FIXED]);
		fscanf(file, "%lg ", &params->system[i].spin[1].magnitude);
		fscanf(file, "%lg ", &params->system[i].spin[1].inclination[FIXED]);
		params->system[i].spin[1].inclination[FIXED] = radianFromDegree(
			params->system[i].spin[1].inclination[FIXED]);
		fscanf(file, "%lg ", &params->system[i].spin[1].azimuth[FIXED]);
		params->system[i].spin[1].azimuth[FIXED] = radianFromDegree(
			params->system[i].spin[1].azimuth[FIXED]);
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
	memset(params, 0, sizeof(SystemParameter));
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
	sscanf(line, "%25s %25s %*s\n", params->amplitude[0], params->amplitude[1]);
	fgets(line, length, file);
	//read_Binary_Parameter_Limits(file, params->system);
}

void readProgramParameters(FILE *file, ProgramParameter *params) {
	assert(file);
	assert(params);
	memset(params, 0, sizeof(ProgramParameter));
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

void printProgramParameters(FILE *file, ProgramParameter *params) {
	fprintf(file, "%10s %10ld\n", "numOfRuns", params->numberOfRuns);
	fprintf(file, "%10s %10hd\n", "prec", params->precision[TO_PLOT]);
	fprintf(file, "%10s %10d\n", "width", params->width[TO_PLOT]);
	fprintf(file, "%10s %10hd\n", "precPlot", params->precision[TO_BACKUP]);
	fprintf(file, "%10s %10d\n", "widthPlot", params->width[TO_BACKUP]);
	fprintf(file, "%10s %10s\n", "folder", params->outputDirectory);
}

void printSystemParameters(FILE *file, SystemParameter *params, OutputFormat *format) {
	fprintf(file, "%10s %10.4lg\n", "freq_I", params->initialFrequency);
	fprintf(file, "%10s %10.4lg\n", "freq_S", params->samplingFrequency);
	fprintf(file, "%10s %10.4lg\n", "time_S", params->samplingTime);
	for (ushort i = 0; i < NUMBER_OF_SYSTEMS; i++) {
		fprintf(file, "%10s %10s\n", "name", params->name[i]);
		fprintf(file, "%10s %10s\n", "approx", params->approximant[i]);
		fprintf(file, "%10s %10s\n", "phase", params->phase[i]);
		fprintf(file, "%10s %10s\n", "spin", params->spin[i]);
		fprintf(file, "%10s %10s\n", "amp", params->amplitude[0]);
		printBinarySystemParameters(file, &params->system[i], format);
	}
}
