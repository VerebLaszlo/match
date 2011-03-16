/**
 * @file util.h
 * @author László Veréb
 * @date 2010.03.27
 * @brief Contains some useful functions.
 */

#ifndef UTIL_H
#define UTIL_H

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#ifndef NDEBUG
#include <string.h>
#endif

typedef enum OUTPUT_CONSTANTS {
	SPECIAL_CHARACTER_LENGTH = 6, SEPARATOR_LENGTH = 3, FORMAT_LENGTH = 10,
} OUTPUT_CONSTANTS;

typedef struct OUTPUT_FORMAT_CONSTANTS {
	short precision;
	short width_Of_Number;
	short width_Of_Number_Width_Separator;
	char separate_Character;
	short left_Justify;
	char format_For_One_Number[FORMAT_LENGTH];
	short precision_To_Plot;
	short width_Of_Number_To_Plot;
	short width_Of_Number_To_Plot_Width_Separator;
	char separate_Character_To_Plot;
	short left_Justify_To_Plot;
	char format_For_One_Number_To_Plot[FORMAT_LENGTH];
	char output_Folder[FILENAME_MAX];
} OUTPUT_FORMAT_CONSTANTS;

/**		Opens the file with the given access mode. On error it terminates the program and prints an
 * error message.
 * @param[in]	file_Name	: the file's path name relative to the program
 * @param[in]	mode		: the mode with which the file is opened
 * @return	pointer to the opened file
 */
FILE * safely_Open_File(char *file_Name, char * mode);

/**		Opens the file just for reading. On error it terminates the program and prints an error
 * message.
 * @param[in]	file_Name	: the file's path name relative to the program
 * @return	pointer to the opened file
 */
FILE * safely_Open_File_For_Reading(char *file_Name);

/**		Opens the file just for writing. On error it terminates the program and prints an error
 * message.
 * @param[in]	file_Name	: the file's path name relative to the program
 * @return	pointer to the opened file
 */
FILE * safely_Open_File_For_Writing(char *file_Name);

/**	Sets the format string for specific number of output.
 * @param[out] format
 * @param[in] number
 * @param[in] format_Constants
 */
void set_Format_For(char format[], const unsigned short number,
		OUTPUT_FORMAT_CONSTANTS*format_Constants);

/**	Sets the format string for specific number of output to plotting.
 * @param[out] format
 * @param[in] number
 * @param[in] format_Constants
 */
void set_Plot_Format_For(char format[], const unsigned short number,
		OUTPUT_FORMAT_CONSTANTS*format_Constants);

/**	Sets the format parameters.
 * @param[in,out] format
 * @param[in] precision
 * @param[in] separator
 * @param[in] left_Justify
 * @param[in] precision_To_Plot
 * @param[in] separator_To_Plot
 * @param[in] left_Justify_To_Plot
 */
void set_Format_Parameters(OUTPUT_FORMAT_CONSTANTS *format, const unsigned short precision,
		const char separator, const unsigned short left_Justify,
		const unsigned short precision_To_Plot, const char separator_To_Plot,
		const unsigned short left_Justify_To_Plot);

/**	Sets the format according to the parameters.
 * @param[in,out] format
 */
void set_Format_For_One_Number(OUTPUT_FORMAT_CONSTANTS *format);

#endif // UTIL_H
