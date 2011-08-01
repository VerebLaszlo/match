/**
 * @file util.h
 * @author László Veréb
 * @date 2011.07.19.
 * @brief Contains useful functions.
 */

#ifndef UTIL_H
#define UTIL_H

#include <stdbool.h>

typedef char *string;	///< shorthand for dynamic string type
typedef unsigned short ushort;	///< shorthand for unsigned short int type

/**	Negates the boolean variable.
 * @param[in,out] var	: boolean variable to be negated.
 */
void neg(bool *var);

#ifdef TEST

/*
typedef enum {
	OFF, BOLD, UNDERSCORE = 4, BLINK, REVERSE, CONCEALED,
} attributes;

typedef enum {
	BLACK = 30, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE,
} foreground;

typedef enum {
	BG_BLACK = 40, BG_RED, BG_GREEN, BG_YELLOW, BG_BLUE, BG_MAGENTA, BG_CYAN, BG_WHITE,
} background;
*/

extern char previous_function[];	///< contains the called function's name
extern char function_file[];		///< contains the called function's filename
extern ushort function_line;		///< contains the calling line of the function
extern char normal[];	///< color code for normal messages
extern char err[];	///< color code for error messages
extern char errBold[];	///< bold color code for error messages
extern char ok[];	///< color code for ok messages
extern char okBold[];	///< bold color code for ok messages

/** Saves the function name to use afterwards. Use it to know the called function's name in the
 * caller. */
#define SET_FUNCTION_NAME() \
	strcpy(previous_function, __func__);

/** Saves the file name to use afterwards. Use it to know the called function's filename in the
 * caller. */
#define SET_FUNCTION_FILE() \
		strcpy(function_file, __FILE__);

/** Saves the file and function name to use afterwards. Use it to know the called function's name
 * and filename in the caller. */
#define SET_FUNCTION_FILE_AND_NAME() \
		SET_FUNCTION_NAME(); \
		SET_FUNCTION_FILE()

/** Saves the line number to use afterwards. Use it before function call, to know where it is
 * called. */
#define SET_FUNCTION_LINE() \
	function_line = __LINE__ + 1;

/**	Prints where was the error.
 */
#define PRINT_ERROR() \
	fprintf(stderr, "%sThe \"%s%s()%s\" function has ERROR!!!\n", err, errBold, previous_function, err); \
	fprintf(stderr, "The error was detected in \"%s%s%s\" at %s%d%s line in \"%s%s()%s\" function from \"%s%s()%s\" function in \"%s%s%s\" file at \"%s%d%s\" line.%s\n", \
		errBold, __FILE__, err, errBold, __LINE__, err, errBold, __func__, err, errBold, previous_function, err, errBold, function_file, err, errBold, function_line, err , normal)

/**	Prints where was the recursive error.
 */
#define PRINT_ERROR_RECURSIVE() \
	fprintf(stderr, "%sRecursive error was in \"%s%s%s\" at %s%d%s line in \"%s%s()%s\" function "\
			"from \"%s%s()%s\".%s\n", err, errBold, __FILE__, err, errBold, __LINE__ - 1, err, errBold, \
			__func__, err, errBold, previous_function, err, normal)

/**	Prints that everything is OK.
 */
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
