/**
 * @file test.h
 * @author László Veréb
 * @date 2011.08.07.
 * @brief Macros for testing functions.
 */

#ifndef TEST_H_
#define TEST_H_

#ifdef TEST
#include <util.h>
#include <string.h>
#include <stdio.h>

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

/**	Structure to contain informations about the called and caller function.
 */
typedef struct {
	char function[FILENAME_MAX];	///< name of the called function
	char file[FILENAME_MAX];///< file of the called function
	ushort line;///< line of the definition of the called function
	char callerFunction[FILENAME_MAX];///< name of the caller function
	char callerFile[FILENAME_MAX];///< file of the caller function
	ushort callerLine;///< line of the calling
}_sourceLocation;

extern _sourceLocation sourceLocation;///< Contains information about the actually tested function.
extern char previous_function[];///< contains the called function's name
extern char function_file[];///< contains the called function's filename
extern ushort function_line;///< contains the calling line of the function
extern char normal[];///< color code for normal messages
extern char err[];///< color code for error messages
extern char errBold[];///< bold color code for error messages
extern char ok[];///< color code for ok messages
extern char okBold[];///< bold color code for ok messages

#define BACKUP_DEFINITION_LINE()\
	static size_t __line_of_definition__;\
	__line_of_definition__ = __LINE__ - 2;

/** Saves the function's name, the name of the file containing it and the line number of the
 * definition of the function. Put it right after the declaration, before everything else!!! */
#define SAVE_FUNCTION_FOR_TESTING() \
	sourceLocation.line = __line_of_definition__;\
	strcpy(sourceLocation.function, __func__);\
	strcpy(sourceLocation.file, __FILE__);

/** Saves the caller function name, the file containing the call and the line where it was called.
 * Put it right before the call!!! */
#define SAVE_FUNCTION_CALLER() \
	strcpy(sourceLocation.callerFunction, __func__);\
	strcpy(sourceLocation.callerFile, __FILE__);\
	sourceLocation.callerLine = __LINE__ + 1;

/**	Prints where was the error.
 */
#define PRINT_ERROR() \
	fprintf(stderr, "%sThe \"%s%s()%s\" function in \"%s%s%s\" at %s%d%s line has %sERROR%s!!!\n",\
			err, errBold, sourceLocation.function, err, errBold, sourceLocation.file, err, errBold,\
			sourceLocation.line, err, errBold, err); \
	fprintf(stderr, "The function was called in \"%s%s%s\" at %s%d%s line from the \"%s%s()%s\""\
			"function,\ndetected at line %s%d%s.%s\n", errBold, sourceLocation.callerFile, err, errBold,\
			sourceLocation.callerLine, err, errBold, sourceLocation.callerFunction, err, errBold, __LINE__, err, normal);

/**	Prints that some of the function are wrong.
 */
#define PRINT_ERROR_FILE() \
	fprintf(stderr, "%sSome of the functions in the \"%s%s%s\" file have %sERRORS%s.%s\n", err,\
			errBold, __FILE__, err, errBold, err, normal)

/**	Prints that everything is OK.
 */
#define PRINT_OK() \
	fprintf(stderr, "%sThe \"%s%s()%s\" function is %sOK%s.%s\n", ok, okBold,\
			sourceLocation.function, ok, okBold, ok, normal)

/**	Prints that everything is OK.
 */
#define PRINT_OK_FILE() \
	fprintf(stderr, "%sEverything in the \"%s%s%s\" file is %sOK%s.%s\n", ok, okBold, __FILE__, ok,\
			okBold, ok, normal)

#else

#define BACKUP_DEFINITION_LINE()
#define SAVE_FUNCTION_FOR_TESTING()
#define SAVE_FUNCTION_CALLER()
#define PRINT_ERROR()
#define PRINT_ERROR_FILE()
#define PRINT_OK()
#define PRINT_OK_FILE()

#endif

#endif /* TEST_H_ */
