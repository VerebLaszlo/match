/**
 * @file util.c
 * @author László Veréb
 * @date 2011.07.19.
 * @brief Contains useful functions.
 */

#include <errno.h>
#include <assert.h>
#include <stdlib.h>
#ifndef NDEBUG
#include <string.h>
#endif
#include "util.h"

extern char * program_invocation_short_name;
extern char * program_invocation_name;

FILE * safelyOpenFile(const char *fileName, const char *mode) {
	assert(strcmp(fileName, ""));
	assert(strcmp(mode, ""));
	FILE *stream;
	errno = 0;
	stream = fopen(fileName, mode);
	if (stream == NULL) {
		if (strcmp(mode, "r")) {
			fprintf(stderr, "%s: Couldn't open file %s for reading; %s\n",
					program_invocation_short_name, fileName, strerror(errno));
		} else if (strcmp(mode, "w")) {
			fprintf(stderr, "%s: Couldn't open file %s for writing; %s\n",
					program_invocation_short_name, fileName, strerror(errno));
		} else if (strcmp(mode, "a")) {
			fprintf(stderr, "%s: Couldn't open file %s for append; %s\n",
					program_invocation_short_name, fileName, strerror(errno));
		} else {
			fprintf(stderr, "%s: Couldn't open file %s; %s\n", program_invocation_short_name,
					fileName, strerror(errno));
		}
		exit(EXIT_FAILURE);
	} else {
		return stream;
	}
}

FILE *safelyOpenForReading(const char *fileName) {
	assert(strcmp(fileName, ""));
	return safelyOpenFile(fileName, "r");
}

FILE *safelyOpenForWriting(const char *fileName) {
	assert(strcmp(fileName, ""));
	return safelyOpenFile(fileName, "w");
}

FILE *safelyOpenForAppend(const char *fileName) {
	assert(strcmp(fileName, ""));
	return safelyOpenFile(fileName, "a");
}

/**	Sets the format string for one number.
 * @param[in,ou]	format	: the format
 */
static void setFormatForOneNumber(OutputFormat *format) {
	assert(format);
	if (format->leftJustified) {
		sprintf(format->oneNumber, "%%- %d.%dlg", format->width, format->precision);
	} else {
		sprintf(format->oneNumber, "%% %d.%dlg", format->width, format->precision);
	}
}

/**	Sets the format variables.
 * @param[out]	format			: format specific variables
 * @param[in]	precision		: user supplied precision of the format
 * @param[in]	width			: user supplied width of the format, its has a minimal value
 * @param[in]	separator		: user defined separator character
 * @param[in]	leftJustified	: left (true) or right (false) justified text
 * @param[in]	name			: user defined name of the format
 * @param[in]	code			: user defined code of the format
 */
static void setOutputFormat(OutputFormat *format, const ushort precision, const ushort width,
		const char separator, bool leftJustified, nameString name, const ushort code) {
	assert(format);
	assert(precision < 100);
	assert(separator);
	assert(name);
	format->precision = precision;
	format->width = width > format->precision + SPECIAL_CHARACTER_LENGTH ? width
			: format->precision + SPECIAL_CHARACTER_LENGTH;
	format->widthWithSeparator = format->width + SEPARATOR_LENGTH;
	format->separator = separator;
	format->leftJustified = leftJustified;
	strcpy(format->name, name);
	format->code = code;
	setFormatForOneNumber(format);
}

/**	Set the format for given number of floating pint data to display.
 * @param[out] formatString	: the generated format string
 * @param[in] number		: number of the data
 * @param[in] format		: contains the format variables
 */
static void setFormat(char formatString[], const ushort number, OutputFormat *format) {
	assert(formatString);
	assert(number);
	assert(format);
	char temp[number * format->widthWithSeparator];
	strcpy(formatString, format->oneNumber);
	for (ushort i = 1; i < number; i++) {
		sprintf(temp, "%s %%%c %s", formatString, format->separator, format->oneNumber);
		strcpy(formatString, temp);
	}
}

/// @name OLD
///@{

void set_Format_For(char format[], const unsigned short number,
		OUTPUT_FORMAT_CONSTANTS*format_Constants) {
	assert(number);
	char temp[number * format_Constants->width_Of_Number_Width_Separator];
	sprintf(format, "%s", format_Constants->format_For_One_Number);
	if (number == 1) {
		return;
	}
	for (short i = 1; i < number; i++) {
		sprintf(temp, "%s %%%c %s", format, format_Constants->separate_Character,
				format_Constants->format_For_One_Number);
		sprintf(format, "%s", temp);
	}
}

void set_Plot_Format_For(char format[], const unsigned short number,
		OUTPUT_FORMAT_CONSTANTS*format_Constants) {
	assert(number);
	char temp[number * format_Constants->width_Of_Number_To_Plot_Width_Separator];
	sprintf(format, "%s", format_Constants->format_For_One_Number_To_Plot);
	if (number == 1) {
		return;
	}
	for (short i = 1; i < number; i++) {
		sprintf(temp, "%s %%%c %s", format, format_Constants->separate_Character_To_Plot,
				format_Constants->format_For_One_Number_To_Plot);
		sprintf(format, "%s", temp);
	}
}

void set_Format_Parameters(OUTPUT_FORMAT_CONSTANTS *format, const unsigned short precision,
		const char separator, const unsigned short left_Justify,
		const unsigned short precision_To_Plot, const char separator_To_Plot,
		const unsigned short left_Justify_To_Plot) {
	assert(precision <100);
	assert(precision_To_Plot < 100);
	format->precision = precision;
	format->width_Of_Number = format->precision + SPECIAL_CHARACTER_LENGTH;
	format->width_Of_Number_Width_Separator = format->width_Of_Number + SEPARATOR_LENGTH;
	format->separate_Character = separator;
	format->left_Justify = left_Justify;
	format->precision_To_Plot = precision_To_Plot;
	format->width_Of_Number_To_Plot = format->precision_To_Plot + SPECIAL_CHARACTER_LENGTH;
	format->width_Of_Number_To_Plot_Width_Separator = format->width_Of_Number_To_Plot
			+ SEPARATOR_LENGTH;
	format->separate_Character_To_Plot = separator_To_Plot;
	format->left_Justify_To_Plot = left_Justify_To_Plot;
	set_Format_For_One_Number(format);
}

inline void set_Format_For_One_Number(OUTPUT_FORMAT_CONSTANTS *format) {
	if (format->left_Justify) {
		sprintf(format->format_For_One_Number, "%%- %d.%dlg", format->width_Of_Number,
				format->precision);
	} else {
		sprintf(format->format_For_One_Number, "%% %d.%dlg", format->width_Of_Number,
				format->precision);
	}
	if (format->left_Justify_To_Plot) {
		sprintf(format->format_For_One_Number_To_Plot, "%%- %d.%dlg",
				format->width_Of_Number_To_Plot, format->precision_To_Plot);
	} else {
		sprintf(format->format_For_One_Number_To_Plot, "%% %d.%dlg",
				format->width_Of_Number_To_Plot, format->precision_To_Plot);
	}
}

///@}
