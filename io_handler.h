/**
 * @file io_handler.h
 *
 * @date Mar 2, 2011
 * @author vereb
 */

#ifndef IO_HANDLER_H_
#define IO_HANDLER_H_

#include "match.h"
#include "generator.h"

extern short is_First;///<a
extern long db;

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
	short amp_Code;
	char approx[2][FILE_NAME_LENGTH];///<a
	char phase[2][FILE_NAME_LENGTH];///<a
	char spin[2][FILE_NAME_LENGTH];///<a
} System_Parameters;

/**
 * Done.
 * @param parameters
 * @param params
 * @param file_Name
 */
void read_Program_Parameters(Program_Parameters *parameters, System_Parameters *params,
		char *file_Name);

/**
 * Done.
 * @param parameters
 * @param file_Name
 */
void read_Parameters(binary_System *parameters, char *file_Name);

/**
 * Done
 * @param prog
 * @param parameters
 * @param sig
 */
void write_Waves_To_Files(Program_Parameters *prog, System_Parameters *parameters,
		signalStruct *sig);

void write_Waves(Program_Parameters *prog, System_Parameters *parameters, signalStruct *sig,
		char *file_Name);

/**
 * Done.
 * @todo előbb megnyitni valahol, hogy megkezdjük írással a fájlt, és ne csak hozzáfüzzünk!!!
 * @param prog
 * @param parameters
 * @param sig
 * @param index
 */
void write_Wave_To_File(Program_Parameters *prog, System_Parameters *parameters, signalStruct *sig,
		short index);

/**
 *
 * @param prog
 * @param parameters
 * @param file_Name
 */
void write_Params_To_File(Program_Parameters *prog, System_Parameters *parameters, char *file_Name);

#endif /* IO_HANDLER_H_ */
