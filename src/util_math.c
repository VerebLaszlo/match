/**
 * @file util_math.c
 * @date 2011.07.19.
 * @author László Veréb
 * @brief Various utilities for math.
 */

#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include "util_math.h"

/// @name Random numbers
///@{

double randomBetweenZeroAndOne(void) {
	return (double) rand() / ((double) RAND_MAX + 1.);
}

double randomBetweenZeroAnd(double top) {
	return top * randomBetweenZeroAndOne();
}

double randomBetween(double bottom, double top) {
	return (top - bottom) * randomBetweenZeroAndOne() + bottom;
}

///@}

/// @name OLD
///@{

const TIME_CONVERSION_CONSTANTS TIME_CONVERSION_CONSTANT = { 60.0, 1.0 / 60.0, //
	60.0, 1.0 / 60.0,//
	60.0 * 60.0, 1.0 / (60.0 * 60.0),//
	24.0, 1.0 / 24.0,//
	24.0 * 60.0, 1.0 / (24.0 * 60.0),//
	24.0 * 60.0 * 60.0, 1.0 / (24.0 * 60.0 * 60.0), };

const CONVERSION_CONSTANTS CONVERSION_CONSTANT = { 180.0 / M_PI, M_PI / 180.0,//
	15.0 * M_PI / 180.0, 1.0 / (15.0 * M_PI / 180.0),//
	15.0, 1.0 / 15.0,//
	15.0 * M_PI / 180.0 / 60.0, 1.0 / (15.0 * M_PI / 180.0 / 60.0),//
	15.0 * M_PI / 180.0 / 60.0 / 60.0, 1.0 / (15.0 * M_PI / 180.0 / 60.0 / 60.0), };

inline long greatest_Number_That_Less_Than(double number) {
	assert(number>0.);
	return (long) exp2(ceil(log(number) / M_LN2));
}

inline long least_Number_That_Greater_Than(double number) {
	assert(number>0.);
	return (long) exp2(floor(log(number) / M_LN2));
}

inline double convert_Time_To_Degree(double hour, double minute, double second) {
	return convert_Time_To_Radian(hour, minute, second) * CONVERSION_CONSTANT.RADIAN_TO_DEGREE;
}

inline double convert_Time_To_Radian(double hour, double minute, double second) {
	return (hour * TIME_CONVERSION_CONSTANT.HOUR_TO_SECOND + minute
			* TIME_CONVERSION_CONSTANT.MINUTE_TO_SECOND + second)
			* CONVERSION_CONSTANT.SECOND_TO_RADIAN;
}

// trigonometric functions

const double epsilon = 1.0e-14;

inline double cos_good(double number) {
	number = normalise_Radians(number);
	if (is_Near(number, M_PI_2, epsilon) || is_Near(number, 3.0 * M_PI_2, epsilon)) {
		return 0.0;
	}
	return cos(number);
}

inline double sin_good(double number) {
	number = normalise_Radians(number);
	if (is_Near(number, M_PI, epsilon) || is_Near(number, 2.0 * M_PI, epsilon)) {
		return 0.0;
	}
	return sin(number);
}

inline double tan_good(double number) {
	number = normalise_Radians(number);
	if (is_Near(number, 0.0, epsilon) || is_Near(number, M_PI, epsilon)) {
		return 0.0;
	}
	if (is_Near(number, M_PI_4, epsilon) || is_Near(number, M_PI + M_PI_4, epsilon)) {
		return 1.0;
	}
	if (is_Near(number, M_PI_2 + M_PI_4, epsilon) || is_Near(number, M_PI + M_PI_2 + M_PI_4,
			epsilon)) {
		return -1.0;
	}
	if (is_Near(number, M_PI_2, epsilon) || is_Near(number, M_PI + M_PI_2, epsilon)) {
		return NAN;
	}
	return tan(number);
}

inline double normalise_Radians(double number) {
	while (number < 0.0) {
		number += 2.0 * M_PI;
	}
	while (number >= 2.0 * M_PI) {
		number -= 2.0 * M_PI;
	}
	return number;
}

inline double is_Near(const double first, const double second, const double epsilon) {
	return fabs(first - second) < epsilon;
}

inline double is_Equal(const double first, const double second) {
	return !islessgreater(first, second);
}

inline double is_Not_Equal(const double first, const double second) {
	return islessgreater(first, second);
}

///@}
