/**
 * @file parameters.h
 *
 * @date Apr 9, 2011
 * @author vereb
 */

#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include "binary_system.h"

typedef enum ParameterConstants_ {
	TO_PLOT, TO_BACKUP, NUMBER_OF_FORMATS, NUMBER_OF_SYSTEMS = 2, LENGTH_OF_STRING = 100,
} ParameterConstants;

typedef struct SystemParameter_ {
	BinarySystem system[NUMBER_OF_SYSTEMS];
	DetectorParameters detector[NUMBER_OF_SYSTEMS];
	double samplingFrequency;
	double samplingTime;
	double initialFrequency;
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

#include "generator.h"

/** System aprameters.
 */
typedef struct System_Parameters {
	binary_System system[2]; ///<a
	double freq_Sampling; ///<a
	double freq_Initial; ///<a
	double time_Sampling; ///<a
	double match_Typ; ///<aa
	double match_Best; ///<a
	double match_Minimax; ///<a
	short shorter; ///<a
	long min_Length; ///<a
	long max_Length; ///<a
	double freq_Min; ///<a
	double freq_Step; ///<a
	double critical_Match; ///<a
	short amp_Code[2]; ///< amplitude code
	char approx[2][FILENAME_MAX]; ///<a
	char phase[2][FILENAME_MAX]; ///<a
	char spin[2][FILENAME_MAX]; ///<a
	char name[2][FILENAME_MAX];
} System_Parameters;

/** Program parameters
 */
typedef struct Program_Parameters {
	char (*output_Directories)[FILENAME_MAX]; ///<a
	long number_Of_Runs; ///<a
	short precision; ///<a
	short precision_To_Plot; ///<a
	short width_Of_Number; ///<a
	short width_Of_Number_To_Plot; ///<a
	char folder[FILENAME_MAX]; ///<a
	double max_Spin; ///<a
	double spin_Step; ///<a
	double freq_Max; ///<a
	double delta_Length; ///<a
	double min_Match; ///<a
} ProgramParameters;

/** Reads the exact system parameters.
 * @param file
 * @param params
 */
void readExactParameters(FILE *file, System_Parameters *params);

/** Reads program parameters
 * @param file
 * @param params
 */
void readProgramParameters(FILE*file, ProgramParameters *params);

/** Prints program parameters.
 * @param file
 * @param params
 */
void print_Program_Parameters(FILE*file, ProgramParameters *params);

/** Reads system parameters.
 * @param file
 * @param params
 */
void readSystemParameters(FILE *file, System_Parameters *params);

/** Prints system parameters.
 * @param file
 * @param params
 */
void print_System_Parameters(FILE *file, System_Parameters *params);

/**
 *
 * @param file
 * @param params
 */
void print_System_Parameters_For_Plot(FILE *file, System_Parameters *params);

#endif /* PARAMETERS_H_ */
