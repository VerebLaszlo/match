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

#endif // UTIL_H
