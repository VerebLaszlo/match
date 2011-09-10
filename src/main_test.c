/**
 * @file main_test.c
 *
 * @date Jul 30, 2011
 * @author vereb
 * @brief
 */

#include <time.h>
#include <math.h>
#include "lal_wrapper.h"

int main(int argc, char *argv[]) {
	printf("%d: %s\n", argc, argv[0]);
#ifdef TEST
	srand(86);
	areUtilMathFunctionsOK();
	areIOFunctionsGood();
	areBinarySystemMassFunctionsGood();
	areBinarySystemSpinFunctionsGood();
	areBinarySystemFunctionsGood();
	areDetectorFunctionsGood();
#endif // TEST
	FILE *file;
	SystemParameter system;
	ProgramParameter program;
	file = safelyOpenForReading(argv[1]);
	readProgramParameters(file, &program);
	//printProgramParameters(stdout, &program);
	fclose(file);
	file = safelyOpenForReading(argv[2]);
	readExactParameters(file, &system);
	convertSpinGlobal(system.system[0].spin);
	convertSpinGlobal(system.system[1].spin);
	//printSystemParameters(stdout, &system, defaultFormat);
	fclose(file);
	SignalStruct signal;
	generateWaveformPair(&system, &signal);
	for (size_t i = 0; i < signal.size; i++) {
		signal.inTime[H1][i] = M_SQRT1_2
			* (signal.componentsInTime[H1P][i] + signal.componentsInTime[H1C][i]);
		signal.inTime[H2][i] = M_SQRT1_2
			* (signal.componentsInTime[H2P][i] + signal.componentsInTime[H2C][i]);
	}
	char fileName[1000];
	sprintf(fileName, "%s/%s", program.outputDirectory, "proba.dat");
	file = safelyOpenForWriting(fileName);
	printTwoSignals(file, &signal, defaultFormat);
	puts("\nOK");
	return 0;
}
