/**
 * @file binary_system.c
 *
 * @date Jul 20, 2011
 * @author vereb
 * @brief Binary system specific
 */

#include "binary_system.h"

/// @name Generation functions
///@{

void generateBinarySystemParameters(binarySystem *system, binarySystem limits[],
		generationMode genMass, generationMode genSpin) {
	BACKUP_DEFINITION_LINE(); //
	assert(system);
	assert(limits);
	system->inclination = randomBetween(limits[MIN].inclination, limits[MAX].inclination);
	system->distance = randomBetween(limits[MIN].distance, limits[MAX].distance);
	for (short i = 0; i < NUMBER_OF_BLACKHOLES; i++) {
		system->flatness[i] = randomBetween(limits[MIN].flatness[i], limits[MAX].flatness[i]);
	}
	massParameters mass[2] = { limits[MIN].mass, limits[MAX].mass };
	generateMass(&system->mass, mass, genMass);
	spinParameters spin[2];
	for (short i = 0; i < NUMBER_OF_BLACKHOLES; i++) {
		spin[MIN] = limits[MIN].spin[0];
		spin[MAX] = limits[MAX].spin[0];
		generateSpin(&system->spin[i], spin, system->inclination, genSpin);
	}
	SAVE_FUNCTION_FOR_TESTING();
}

///@}
/// @name Printing functions
///@{

void printBinarySystemParameters(FILE *file, binarySystem *system, OutputFormat *format) {
	BACKUP_DEFINITION_LINE();
	printMassParameters(file, &system->mass, format);
	printSpinParameters(file, &system->spin[0], format);
	printSpinParameters(file, &system->spin[1], format);
	ushort number = 3;
	char formatString[number * format->widthWithSeparator];setFormat
	(formatString, number, format);
	fprintf(file, formatString, system->flatness[0], system->flatness[1], system->inclination);
	setFormatEnd(formatString, number, format);
	fprintf(file, formatString, system->distance, system->coalescencePhase,
			system->coalescenceTime);
	SAVE_FUNCTION_FOR_TESTING();
}

///@}

#ifdef TEST
/// @name Testing functions
///@{

bool areBinarySystemFunctionsGood(void) {
	return true;
}

///@}
#endif	// TEST
