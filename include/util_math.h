/**
 * @file util_math.h
 * @date 2011.07.19.
 * @author László Veréb
 * @brief Various utilities for math.
 */

#ifndef UTIL_MATH_H
#define UTIL_MATH_H

#include "util.h"

/**	Inserts the square of the number.
 * @param[in] number
 */
#define square(square_number) \
	({ typeof(square_number) _square_number = (square_number); _square_number*_square_number; })

/**	Inserts the cube of the number.
 * @param[in] number
 */
#define cube(cube_number) \
	({ typeof(cube_number) _cube_number = (cube_number); square(_cube_number)*_cube_number; })

/// @name Trigonometric functions
///@{

/**	Better sine function. Returns \f$0\f$ at \f$\pi, 2\pi\f$. @warning Use the normaliseRadians()
 * function beforehand, if the argument is not between \f$[0, 2\pi)\f$.
 * @param[in] number : argument
 * @return \f$\sin(number)\f$
 */
double sinGood(double number);

/**	Better cosine function. Returns \f$0\f$ at \f$\pi/2, 3\pi/2\f$. @warning Use the normaliseRadians()
 * function beforehand, if the argument is not between \f$[0, 2\pi)\f$.
 * @param[in] number : argument
 * @return \f$\cos(number)\f$
 */
double cosGood(double number);

/**	Better tangent function. Returns \f$0\f$ at \f$\pi\f$, and \f$\infty\f$ at \f$\pi/2, 3\pi/2\f$.
 * @warning Use the normaliseRadians() function beforehand, if the argument is not between
 * \f$[0, 2\pi)\f$.
 * @param[in] number : argument
 * @return \f$\tan(number)\f$
 */
double tanGood(double number);

/**	Returns the normalized angle between \f$[0,2\pi)\f$.
 * @param[in] number : argument
 * @return normalized angle
 */
double normaliseRadians(double number);

/**	Converts angle units.
 * @param[in] degree
 * @return
 */
double radianFromDegree(double degree);

/**	Converts angle units.
 * @param[in] turn
 * @return
 */
double radianFromTurn(double turn);

/**	Converts angle units.
 * @param[in] radian
 * @return
 */
double turnFromRadian(double radian);

/**	Converts angle units.
 * @param[in] degree
 * @return
 */
double turnFromDegree(double degree);

/**	Converts angle units.
 * @param[in] radian
 * @return
 */
double degreeFromRadian(double radian);

/**	Converts angle units.
 * @param[in] turn
 * @return
 */
double degreeFromTurn(double turn);

///@}
/// @name Random numbers
///@{

/**  Returns a pseudo-random number in the range [0,1). Use srand() beforehand to initialize the
 * random number generator.
 * @return a random number
 */
double randomBetweenZeroAndOne(void);

/**  Returns a pseudo-random number in the range [0,top). Use srand() beforehand to initialize
 * the random number generator.
 * @param[in]	top	: the upper limit of the generated number.
 * @return a random number.
 */
double randomBetweenZeroAnd(double top);

/**  Returns a pseudo-random number in the range [bottom,top). Use srand() beforehand to
 * initialize the random number generator. It works when bottom is greater than top, too.
 * @param[in]	bottom	: the lowest possible number
 * @param[in]	top		: the upper limit of the generated numbers.
 * @return a random number.
 */
double randomBetween(double bottom, double top);

///@}

extern const double EPSILON;

bool isNear(const double first, const double second, const double epsilon);

#ifdef TEST

bool isOK_randomBetween(void);

bool areUtilMathFunctionsOK(void);

#endif	// TEST
/// @name OLD
///@{

/** Time conversion constants
 */
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

/** Degree conversion constants
 */
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

// trigonometric functions

/**	Returns nonzero value, if the difference is less then epsilon.
 * @param[in] first
 * @param[in] second
 * @param[in] epsilon
 * @return
 */
double is_Near(const double first, const double second, const double epsilon);

/**	Returns nonzero value if the numbers are equal.
 * @param[in] first
 * @param[in] second
 * @return
 */
double is_Equal(const double first, const double second);

/**	Returns nonzero value if the numbers are not equal.
 * @param[in] first
 * @param[in] second
 * @return
 */
double is_Not_Equal(const double first, const double second);

///@}

#endif	// UTIL_MATH_H
