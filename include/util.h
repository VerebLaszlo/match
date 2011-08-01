/**
 * @file util.h
 * @author László Veréb
 * @date 2011.07.19.
 * @brief Contains useful functions.
 */

#ifndef UTIL_H
#define UTIL_H

#include <stdbool.h>

typedef char *string;
typedef unsigned short ushort;

/**	Negates the boolean variable.
 * @param[in,ou] var	: boolean variable to be negated.
 */
void neg(bool *var);

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

extern char previous_function[];
extern char function_file[];
extern ushort function_line;
extern char normal[8];
//char err[8];
//char errBold[8];
//char ok[8];
//char okBold[8];

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
	fprintf(stderr, "%sRecursive error was in \"%s%s%s\" at %s%d%s line in \"%s%s()%s\" function "\
			"from \"%s%s()%s\".%s\n", err, errBold, __FILE__, err, errBold, __LINE__ - 1, err, errBold, \
			__func__, err, errBold, previous_function, err, normal)

#define PRINT_OK() \
	fprintf(stderr, "%sThe \"%s%s()%s\" function is OK.%s\n", ok, okBold, previous_function, ok, normal)

#else
#define SET_FUNCTION_NAME()
#define SET_FUNCTION_FILE()
#define SET_FUNCTION_FILE_AND_NAME()
#define SET_FUNCTION_LINE()
#define PRINT_ERROR()
#define PRINT_ERROR_RECURSIVE()
#define PRINT_OK()
#endif

#endif // UTIL_H
