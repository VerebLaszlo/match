/**
 * @file util.h
 * @author László Veréb
 * @date 2011.07.19.
 * @brief Contains useful functions.
 */

#ifndef UTIL_H
#define UTIL_H

#include <stdio.h>
#include <stdbool.h>

typedef char *string;
typedef unsigned short ushort;

#ifdef TEST

extern char err[8];
extern char errBold[8];
extern char ok[8];
extern char okBold[8];

typedef enum {
	OFF, BOLD, UNDERSCORE = 4, BLINK, REVERSE, CONCEALED,
} attributes;

typedef enum {
	BLACK = 30, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE,
} foreground;

typedef enum {
	BG_BLACK = 40, BG_RED, BG_GREEN, BG_YELLOW, BG_BLUE, BG_MAGENTA, BG_CYAN, BG_WHITE,
} background;

extern char previous_function[FILENAME_MAX];
extern char function_file[FILENAME_MAX];
extern ushort function_line;

#define SET_FUNCTION_NAME() \
	strcpy(previous_function, __func__);

#define SET_FUNCTION_FILE() \
		strcpy(function_file, __FILE__);

#define SET_FUNCTION_FILE_AND_NAME() \
		SET_FUNCTION_NAME(); \
		SET_FUNCTION_FILE()

#define SET_FUNCTION_LINE() \
	function_line = __LINE__ + 1;

#define PRINT_ERROR() \
	fprintf(stderr, "%sThe \"%s%s()%s\" function has ERROR!!!\n", err, errBold, previous_function, err); \
	fprintf(stderr, "The error was detected in \"%s%s%s\" at %s%d%s line in \"%s%s()%s\" function from \"%s%s()%s\" function in \"%s%s%s\" file at \"%s%d%s\" line.%s\n", \
		errBold, __FILE__, err, errBold, __LINE__, err, errBold, __func__, err, errBold, previous_function, err, errBold, function_file, err, errBold, function_line, err , normal)

#define PRINT_ERROR_RECURSIVE() \
	fprintf(stderr, "Recursive error was in \"%s%s%s\" at %s%d%s line in \"%s%s()%s\" function "\
			"from \"%s%s()%s\".%s\n", errBold, __FILE__, err, errBold, __LINE__ - 1, err, errBold, \
			__func__, err, errBold, previous_function, err, normal)

#define PRINT_OK() \
	fprintf(stderr, "%sThe \"%s%s()%s\" function is OK.%s\n", ok, okBold, previous_function, ok, normal)

#else
#define SET_FUNCTION_NAME()
#define SET_FUNCTION_LINE()
#endif

/**	Negates the boolean variable.
 * @param[in,ou] var	: boolean variable to be negated.
 */
void neg(bool *var);

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
/// @name OLD
///@{

/** Specifies the formatting of the output
 */
typedef struct OUTPUT_FORMAT_CONSTANTS {
	ushort precision; ///< precision of the number
	ushort width_Of_Number; ///< length of the string containing the number without separator character
	ushort width_Of_Number_Width_Separator; ///< the length of the string containing the number with separator character
	char separate_Character; ///< the separator character
	short left_Justify; ///< true if the numbers are left justified
	char format_For_One_Number[FORMAT_LENGTH]; ///< contains the format-string for one number
	ushort precision_To_Plot; ///< precision of the number for plotting
	ushort width_Of_Number_To_Plot; ///< length of the string containing the number without separator character for plotting
	ushort width_Of_Number_To_Plot_Width_Separator; ///< the length of the string containing the number with separator character for plotting
	char separate_Character_To_Plot; ///< the separator character for plotting
	ushort left_Justify_To_Plot; ///< true if the numbers are left justified for plotting
	char format_For_One_Number_To_Plot[FORMAT_LENGTH]; ///< contains the format-string for one number for plotting
	char output_Folder[FILENAME_MAX]; ///< where to write files
} OUTPUT_FORMAT_CONSTANTS;

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

///@}

#endif // UTIL_H
