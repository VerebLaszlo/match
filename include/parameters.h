/**
 * @file parameters.h
 *
 * @date Apr 9, 2011
 * @author vereb
 */

#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include "binary_system.h"
#include "detector.h"

typedef enum ParameterConstants_ {
	TO_PLOT, TO_BACKUP, NUMBER_OF_FORMATS, NUMBER_OF_SYSTEMS = 2, LENGTH_OF_STRING = 100,
} ParameterConstants;

typedef struct SystemParameter_ {
	BinarySystem system[NUMBER_OF_SYSTEMS];
	DetectorParameters detector[NUMBER_OF_SYSTEMS];
	double coalescencePhase[NUMBER_OF_SYSTEMS];
	double coalescenceTime[NUMBER_OF_SYSTEMS];
	double samplingFrequency;
	double samplingTime;
	double initialFrequency;
	double endingFrequency;
	char name[NUMBER_OF_SYSTEMS][LENGTH_OF_STRING];
	char approximant[NUMBER_OF_SYSTEMS][LENGTH_OF_STRING];
	char phase[NUMBER_OF_SYSTEMS][LENGTH_OF_STRING];
	char spin[NUMBER_OF_SYSTEMS][LENGTH_OF_STRING];
	char amplitude[NUMBER_OF_SYSTEMS][LENGTH_OF_STRING];
} SystemParameter;

typedef struct ProgramParameter_ {
	char outputDirectory[FILENAME_MAX];
	ulong numberOfRuns;
	ushort precision[NUMBER_OF_FORMATS];
	ushort width[NUMBER_OF_FORMATS];
} ProgramParameter;

void readExactParameters(FILE *file, SystemParameter *params);

void readSystemParameters(FILE *file, SystemParameter *params);

void readProgramParameters(FILE *file, ProgramParameter *params);

void printProgramParameters(FILE *file, ProgramParameter *params);

void printSystemParameters(FILE *file, SystemParameter *params, OutputFormat *format);

#endif /* PARAMETERS_H_ */
