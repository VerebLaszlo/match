/**
 * @file program_functions.c
 * @author vereb
 * @date Sep 27, 2011
 * @brief 
 */

#include <math.h>
#include "program_functions.h"
#include "match.h"
#include "lal_wrapper.h"

static void initializeGeneration(ProgramParameter *program, SystemParameter *parameters,
	char *programFile, char *parametersFile) {
	FILE *file = safelyOpenForReading(programFile);
	readProgramParameters(file, program);
	fclose(file);
	file = safelyOpenForReading(parametersFile);
	readExactParameters(file, parameters);
	convertSpinGlobal(parameters->system[0].spin);
	convertSpinGlobal(parameters->system[1].spin);
	fclose(file);

}

void run(char *programFile, char *parametersFile, bool plot, bool calculateMatches) {
	SystemParameter parameters;
	ProgramParameter program;
	initializeGeneration(&program, &parameters, programFile, parametersFile);
	setSignalExistanceFunctions(calculateMatches);
	SignalStruct signal;
	generateWaveformPair(&parameters, &signal, calculateMatches);
	if (plot) {
		for (size_t i = 0; i < signal.size; i++) {
			signal.inTime[H1][i] = M_SQRT1_2
				* (signal.componentsInTime[H1P][i] + signal.componentsInTime[H1C][i]);
			signal.inTime[H2][i] = M_SQRT1_2
				* (signal.componentsInTime[H2P][i] + signal.componentsInTime[H2C][i]);
		}
		char fileName[1000];
		sprintf(fileName, "%s/%s", program.outputDirectory, "proba.dat");
		FILE *file = safelyOpenForWriting(fileName);
		printTwoSignals(file, &signal, defaultFormat);
		fclose(file);
	}
	if (calculateMatches) {
		double type, minimax, best;
		size_t min, max;
		calculateIndexBoundariesFromFrequencies(parameters.initialFrequency,
			parameters.endingFrequency, parameters.samplingFrequency, &min, &max);
		calc_Matches(&signal, min, max, &type, &best, &minimax);
		printf("%lg %lg %lg\n", minimax, type, best);
	}
	destroySignal(&signal);
}
