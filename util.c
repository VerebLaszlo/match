/*
 * @file util.c
 * @author László Veréb
 * @date 2010.01.11.
 */

#include "util.h"

extern char * program_invocation_short_name;
extern char * program_invocation_name;

FILE * safely_Open_File(char *file_Name, char *mode) {
	assert(strcmp(file_Name, ""));
	assert(strcmp(mode, ""));
	FILE *stream;
	errno = 0;
	stream = fopen(file_Name, mode);
	if (stream == NULL) {
		fprintf(stderr, "%s: Couldn't open file %s; %s\n", program_invocation_short_name,
				file_Name, strerror(errno));
		exit(EXIT_FAILURE);
	} else {
		return stream;
	}
}

FILE * safely_Open_File_For_Reading(char *file_Name) {
	assert(strcmp(file_Name, ""));
	FILE *stream;
	errno = 0;
	stream = fopen(file_Name, "r");
	if (stream == NULL) {
		fprintf(stderr, "%s: Couldn't open file for reading %s; %s\n",
				program_invocation_short_name, file_Name, strerror(errno));
		exit(EXIT_FAILURE);
	} else {
		return stream;
	}
}

FILE * safely_Open_File_For_Writing(char *file_Name) {
	assert(strcmp(file_Name, ""));
	FILE *stream;
	errno = 0;
	stream = fopen(file_Name, "w");
	if (stream == NULL) {
		fprintf(stderr, "%s: Couldn't open file for writing %s; %s\n",
				program_invocation_short_name, file_Name, strerror(errno));
		exit(EXIT_FAILURE);
	} else {
		return stream;
	}
}

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
