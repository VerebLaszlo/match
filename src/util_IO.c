/**
 * @file util_IO.c
 *
 * @date Aug 1, 2011
 * @author vereb
 * @brief Handles input-output specific events.
 */

#include <errno.h>
#include <stdlib.h>
#include "util_IO.h"
#ifndef NDEBUG
#include <assert.h>
#include <string.h>
#endif

extern char * program_invocation_short_name;	///< short name of the program
extern char * program_invocation_name;			///< long name of the program

/// @name File handling functions
///@{

FILE *safelyOpenFile(const char *fileName, const char *mode) {
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

///@}
/// @name Output formatting functions and types
///@{

/**	Sets the format string for one number.
 * @param[in,out]	format	: the format
 */
static void setFormatForOneNumber(OutputFormat *format) {
	SAVE_FUNCTION_FOR_TESTING();
	assert(format);
	assert(format->width > 0);
	if (format->leftJustified) {
		sprintf(format->oneNumber, "%%- %d.%dlg", format->width, format->precision);
	} else {
		sprintf(format->oneNumber, "%% %d.%dlg", format->width, format->precision);
	}
}

void setOutputFormat(OutputFormat *format, const ushort precision, const ushort width,
		const char separator, bool leftJustified, nameString name, const ushort code) {
	SAVE_FUNCTION_FOR_TESTING();
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
}

void setFormat(char formatString[], const ushort number, OutputFormat *format) {
	SAVE_FUNCTION_FOR_TESTING();
	assert(formatString);
	assert(number);
	assert(format);
	char temp[number * format->widthWithSeparator];strcpy
	(formatString, format->oneNumber);
	for (ushort i = 1; i < number; i++) {
		sprintf(temp, "%s %%%c %s", formatString, format->separator, format->oneNumber);
		strcpy(formatString, temp);
	}
	//strcpy(temp, formatString);
	//sprintf(formatString, "%s\n", temp);
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

///@}

#ifdef TEST

/// @name test functions
///@{

static bool isOK_setFormatForOneNumber(void) {
	OutputFormat format;
	nameString result[2] = { "%- 10.5lg", "% 10.5lg" };
	format.leftJustified = true;
	format.precision = 5;
	format.width = 10;
	for (short i = 0; i < 2; i++, neg(&format.leftJustified)) {
		SAVE_FUNCTION_CALLER();
		setFormatForOneNumber(&format);
		if (strcmp(result[i], format.oneNumber)) {
			PRINT_ERROR();
			return false;
		}
	}
	format.leftJustified = false;
	for (short i = 0; i < 2; i++, neg(&format.leftJustified)) {
		SAVE_FUNCTION_CALLER();
		setFormatForOneNumber(&format);
		if (!strcmp(result[i], format.oneNumber)) {
			PRINT_ERROR();
			return false;
		}
	}
	PRINT_OK();
	return true;
}

static bool isOK_setOutputFormat(void) {
	if (!isOK_setFormatForOneNumber()) {
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
			SAVE_FUNCTION_CALLER();
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
			if (strcmp(name[k], format.name)) {
				PRINT_ERROR();
				return false;
			}
			width += 2 * SPECIAL_CHARACTER_LENGTH;
		}
	}
	/// @todo a hibás bemenet nincs leellenőrizve!!!
	PRINT_OK();
	return true;
}

static bool isOK_setFormat(void) {
	if (!isOK_setOutputFormat()) {
		return false;
	}
	OutputFormat format;
	char separator = '%';
	ushort code = 1;
	nameString name = "multiformat";
	ushort precision = 5;
	bool leftJusfified = false;
	ushort width = 12;
	nameString output;
	nameString result[2][2] = { { "% 12.5lg", "% 12.5lg %% % 12.5lg" }, //
			{ "%- 12.5lg", "%- 12.5lg %% %- 12.5lg" } };
	for (short i = 0; i < 2; i++, neg(&leftJusfified)) {
		setOutputFormat(&format, precision, width, separator, leftJusfified, name, code);
		SAVE_FUNCTION_CALLER();
		setFormat(output, 1, &format);
		if (strcmp(result[i][0], output)) {
			PRINT_ERROR();
			return false;
		}
		SAVE_FUNCTION_CALLER();
		setFormat(output, 2, &format);
		if (strcmp(result[i][1], output)) {
			PRINT_ERROR();
			return false;
		}
	}
	PRINT_OK();
	return true;
}

bool areIOFunctionsGood(void) {
	if (isOK_setFormat()) {
		PRINT_OK_FILE();
		return true;
	}
	PRINT_ERROR_FILE();
	return false;
}

///@}

#endif // TEST
