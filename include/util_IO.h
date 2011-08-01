/**
 * @file util_IO.h
 *
 * @date Aug 1, 2011
 * @author vereb
 * @brief
 */

#ifndef UTIL_IO_H_
#define UTIL_IO_H_

#include <stdbool.h>
#include <stdio.h>
#include "util.h"

/// @name File handling functions
///@{

/**	Opens the file with the given access mode. On error it terminates the program and prints an
 * error message.
 * @param[in]	fileName	: the file's path name relative to the program
 * @param[in]	mode		: the mode with which the file is opened
 * @return	pointer to the opened file
 */
FILE *safelyOpenFile(const char *fileName, const char *mode);

/**	Opens the file just for reading. On error it terminates the program and prints an error
 * message.
 * @param[in]	fileName	: the file's path name relative to the program
 * @return	pointer to the opened file
 */
FILE *safelyOpenForReading(const char *file_Name);

/**	Opens the file just for writing. On error it terminates the program and prints an error
 * message.
 * @param[in]	fileName	: the file's path name relative to the program
 * @return	pointer to the opened file
 */
FILE *safelyOpenForWriting(const char *file_Name);

/**	Opens the file just for append. On error it terminates the program and prints an error
 * message.
 * @param[in]	fileName	: the file's path name relative to the program
 * @return	pointer to the opened file
 */
FILE *safelyOpenForAppend(const char *fileName);

///@}
/// @name Output formatting functions and types

/** Constants for output
 */
typedef enum OUTPUT_CONSTANTS {
	SPECIAL_CHARACTER_LENGTH = 6, MAXIMUM_WIDTH = 99, MAXIMUM_PRECISION = MAXIMUM_WIDTH
			- SPECIAL_CHARACTER_LENGTH, SEPARATOR_LENGTH = 3, FORMAT_LENGTH = 11, NAMES_LENGTH = 100,
} OUTPUT_CONSTANTS;

typedef char nameString[NAMES_LENGTH];

/**	Contains values to format an output.
 */
typedef struct tagOutputFormat {
	ushort precision;
	ushort width;
	ushort widthWithSeparator;
	char separator;
	bool leftJustified;
	char oneNumber[FORMAT_LENGTH];
	nameString name;
	ushort code;
} OutputFormat;

/**	Sets the format variables.
 * @param[out]	format			: format specific variables
 * @param[in]	precision		: user supplied precision of the format
 * @param[in]	width			: user supplied width of the format, its has a minimal value
 * @param[in]	separator		: user defined separator character
 * @param[in]	leftJustified	: left (true) or right (false) justified text
 * @param[in]	name			: user defined name of the format
 * @param[in]	code			: user defined code of the format
 */
void setOutputFormat(OutputFormat *format, const ushort precision, const ushort width,
		const char separator, bool leftJustified, nameString name, const ushort code);

/**	Set the format for given number of floating pint data to display.
 * @param[out] formatString	: the generated format string
 * @param[in] number		: number of the data
 * @param[in] format		: contains the format variables
 */
void setFormat(char formatString[], const ushort number, OutputFormat *format);

#ifdef TEST

/// @name test functions
///@{

bool isOK_setFormatForOneNumber(void);

bool isOK_setOutputFormat(void);

///@}

#endif // TEST

///@}

#endif /* UTIL_IO_H_ */
