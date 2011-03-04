/**
 * @file util_math.h
 * @date 2010-08-16 10.27.00
 * @author László Veréb
 * @brief Various utilities for math. @callgraph @callergraph
 */

#ifndef UTIL_MATH_H
#define UTIL_MATH_H

#include <math.h>
#include "util.h"

/**	Calculates the square of the parameter.
 * The argument evaluates only once, so there is no undefined side effect
 * @todo examine the side effects
 * @param[in] number
 */
#define SQR(number) \
	({ typeof(number) _number = (number); _number*_number; })

typedef struct {
	const double MINUTE_TO_SECOND;
	const double SECOND_TO_MINUTE;
	const double HOUR_TO_MINUTE;
	const double MINUTE_TO_HOUR;
	const double HOUR_TO_SECOND;
	const double SECOND_TO_HOUR;
	const double DAY_TO_HOUR;
	const double HOUR_TO_DAY;
	const double DAY_TO_MINUTE;
	const double MINUTE_TO_DAY;
	const double DAY_TO_SECOND;
	const double SECOND_TO_DAY;
} TIME_CONVERSION_CONSTANTS;

typedef struct {
	const double RADIAN_TO_DEGREE;
	const double DEGREE_TO_RADIAN;
	const double HOUR_TO_RADIAN;
	const double RADIAN_TO_HOUR;
	const double HOUR_TO_DEGREE;
	const double DEGREE_TO_HOUR;
	const double MINUTE_TO_RADIAN;
	const double RADIAN_TO_MINUTE;
	const double SECOND_TO_RADIAN;
	const double RADIAN_TO_SECOND;
} CONVERSION_CONSTANTS;

/** Contains constants to convert from and to various time measurement units.
 */
extern const TIME_CONVERSION_CONSTANTS TIME_CONVERSION_CONSTANT;

/** Contains constants to convert from and to various time and degree measurement units.
 */
extern const CONVERSION_CONSTANTS CONVERSION_CONSTANT;

/**  Returns a pseudo-random number in the range [0,top). Use srand() beforehand to initialize
 * the random number generator.
 * @param[in]	top	: the upper limit of the generated number.
 * @return a random number.
 */
double random_From_Zero_To(double top);

/**  Returns a pseudo-random number in the range [bottom,top). Use srand() beforehand to
 * initialize the random number generator. It works when bottom is greater than top, too.
 * @param[in]	bottom	: the lowest possible number
 * @param[in]	top		: the upper limit of the generated numbers.
 * @return a random number.
 */
double random_Between(double bottom, double top);

/**  Returns a pseudo-random number in the range [0,1). Use srand() beforehand to initialize the
 * random number generator.
 * @return a random number
 */
double random_From_Zero_To_One(void);

/**	Returns the smallest integral value that is not less than number and is the power of two.
 * @param[in]	number	: the number in question
 * @return
 */
long greatest_Number_That_Less_Than(double number);

/**	Returns the smallest integral value that is not less than number and is the power of two.
 * @param[in]	number	: the number in question
 * @return
 */
long least_Number_That_Greater_Than(double number);

/**	Converts time into radian.
 * @param[in] hour		: hours
 * @param[in] minute	: minutes
 * @param[in] second	: seconds
 * @return the time value in radian
 */
double convert_Time_To_Radian(double hour, double minute, double second);

/**	Converts time into degree.
 * @param[in] hour		: hours
 * @param[in] minute	: minutes
 * @param[in] second	: seconds
 * @return the time value in degree
 */
double convert_Time_To_Degree(double hour, double minute, double second);

#endif	// UTIL_MATH_H
