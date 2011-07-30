/**
 * @file util.c
 * @author László Veréb
 * @date 2011.07.19.
 * @brief Contains useful functions.
 */

#include <errno.h>
#include <stdlib.h>
#include "util.h"
#ifndef NDEBUG
#include <assert.h>
#include <string.h>
#endif

#ifdef TEST
char previous_function[FILENAME_MAX];
char function_file[FILENAME_MAX];
ushort function_line;
char normal[8] = "\e[0m";
char err[8] = "\e[0;31m";
char errBold[8] = "\e[1;31m";
char ok[8] = "\e[0;36m";
char okBold[8] = "\e[1;36m";
#endif

extern char * program_invocation_short_name;
extern char * program_invocation_name;

void neg(bool *var) {
	*var = !*var;
}

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

static void printFormat(FILE *file, OutputFormat *format) {
	fprintf(file, "prefision:             %u\n", format->precision);
	fprintf(file, "width:                 %u\n", format->width);
	fprintf(file, "width with separator:  %u\n", format->widthWithSeparator);
	fprintf(file, "separator:             %c\n", format->separator);
	fprintf(file, "left jusfified:        %d\n", format->leftJustified);
	fprintf(file, "format for one number: %s\n", format->oneNumber);
	fprintf(file, "name of the format:    %s\n", format->name);
	fprintf(file, "code of the format:    %d\n", format->code);
}

/**	Sets the format string for one number.
 * @param[in,ou]	format	: the format
 */
static void setFormatForOneNumber(OutputFormat *format) {
	assert(format);
	assert(format->width > 0);
	if (format->leftJustified) {
		sprintf(format->oneNumber, "%%- %d.%dlg", format->width, format->precision);
	} else {
		sprintf(format->oneNumber, "%% %d.%dlg", format->width, format->precision);
	}
	SET_FUNCTION_FILE_AND_NAME();
}

void setOutputFormat(OutputFormat *format, const ushort precision, const ushort width,
		const char separator, bool leftJustified, nameString name, const ushort code) {
	assert(format);
	assert(width < MAXIMUM_WIDTH);
	assert(precision < MAXIMUM_PRECISION);
	assert(separator);
	assert(name);
	format->precision = precision;
	format->width =
			width > format->precision + SPECIAL_CHARACTER_LENGTH ? width :
					format->precision + SPECIAL_CHARACTER_LENGTH;
	format->widthWithSeparator = format->width + SEPARATOR_LENGTH;
	format->separator = separator;
	format->leftJustified = leftJustified;
	strcpy(format->name, name);
	format->code = code;
	setFormatForOneNumber(format);
	SET_FUNCTION_FILE_AND_NAME();
}

void setFormat(char formatString[], const ushort number, OutputFormat *format) {
	assert(formatString);
	assert(number);
	assert(format);
	char temp[number * format->widthWithSeparator];strcpy
	(formatString, format->oneNumber);
	for (ushort i = 1; i < number; i++) {
		sprintf(temp, "%s %%%c %s", formatString, format->separator, format->oneNumber);
		strcpy(formatString, temp);
	}
	SET_FUNCTION_FILE_AND_NAME();
}

#ifdef TEST

/// @name test functions
///@{

bool isOK_setFormatForOneNumber(void) {
	OutputFormat format;
	nameString result[2] = { "%- 10.5lg", "% 10.5lg" };
	format.leftJustified = true;
	format.precision = 5;
	format.width = 10;
	for (short i = 0; i < 2; i++, neg(&format.leftJustified)) {
		SET_FUNCTION_LINE();
		setFormatForOneNumber(&format);
		if (strcmp(result[i], format.oneNumber)) {
			PRINT_ERROR();
			return false;
		}
	}
	format.leftJustified = false;
	for (short i = 0; i < 2; i++, neg(&format.leftJustified)) {
		SET_FUNCTION_LINE();
		setFormatForOneNumber(&format);
		if (!strcmp(result[i], format.oneNumber)) {
			PRINT_ERROR();
			return false;
		}
	}
	PRINT_OK();
	SET_FUNCTION_FILE_AND_NAME();
	return true;
}

bool isOK_setOutputFormat(void) {
	if (!isOK_setFormatForOneNumber()) {
		PRINT_ERROR_RECURSIVE();
		return false;
	}
	OutputFormat format;
	char separator = '%';
	ushort code = 0;
	nameString name[6] = { "ab", "bd", "cf", "dh", "ej", "fl" };
	ushort precision = 2 * SPECIAL_CHARACTER_LENGTH;
	bool leftJusfified = false;
	ushort width;
	for (short i = 0; i < 2; i++, neg(&leftJusfified)) {
		width = precision - 2 * SPECIAL_CHARACTER_LENGTH;
		for (short j = 0; j < 3; j++, code++) {
			short k = 3 * i + j;
			SET_FUNCTION_LINE();
			setOutputFormat(&format, precision, width, separator, leftJusfified, name[k], code);
			if (format.precision != precision) {
				PRINT_ERROR();
				return false;
			}
			if (format.width < width) {
				PRINT_ERROR();
				return false;
			}
			if (format.width < precision) {
				PRINT_ERROR();
				return false;
			}
			if (!strcmp(name[k], format.name)) {
				PRINT_ERROR();
				return false;
			}
			printFormat(stdout, &format);
			width += 2 * SPECIAL_CHARACTER_LENGTH;
		}
	}
	PRINT_OK();
	SET_FUNCTION_FILE_AND_NAME();
	return true;
}

///@}

#endif // TEST
/// @name OLD
///@{

void set_Format_For(char format[], const unsigned short number,
		OUTPUT_FORMAT_CONSTANTS*format_Constants) {
	assert(number);
	char temp[number * format_Constants->width_Of_Number_Width_Separator];sprintf
	(format, "%s", format_Constants->format_For_One_Number);
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
	char temp[number * format_Constants->width_Of_Number_To_Plot_Width_Separator];sprintf
	(format, "%s", format_Constants->format_For_One_Number_To_Plot);
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
