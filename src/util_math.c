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

static const double RADIAN_PER_DEGREE = M_PI / 180.0;	///< degree to radian conversion constant
static const double RADIAN_PER_TURN = M_PI + M_PI;	///< turn to radian conversion constant
static const double DEGREE_PER_TURN = 360.0;	///< turn to degree conversion constant

/// @name Trigonometric functions
///@{

double sinGood(double number) {
	if (number == M_PI || number = M_PI + M_PI) {
		return 0.0;
	}
	return sin(number);
}

double cosGood(double number) {
	if (number == M_PI_2 || number = M_PI_2 + M_PI) {
		return 0.0;
	}
	return cos(number);
}

double tanGood(double number) {
	if (number == M_PI) {
		return 0.0;
	}
	if (number == M_PI_2 || number == M_PI_2 + M_PI) {
		return NAN;
	}
	return tan(number);
}

double normaliseRadians(double number) {
	while (number < 0.0) {
		number += M_PI + M_PI;
	}
	while (number >= M_PI + M_PI) {
		number -= M_PI + M_PI;
	}
	return number;
}

double radianFromDegree(double degree) {
	return degree * RADIAN_PER_DEGREE;
}

double radianFromTurn(double turn) {
	return turn * RADIAN_PER_TURN;
}

double turnFromRadian(double radian) {
	return radian / RADIAN_PER_TURN;
}

double turnFromDegree(double degree) {
	return degree / DEGREE_PER_TURN;
}

double degreeFromRadian(double radian) {
	return radian / RADIAN_PER_DEGREE;
}

double degreeFromTurn(double turn) {
	return turn * DEGREE_PER_TURN;
}

///@}
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
60.0, 1.0 / 60.0, //
60.0 * 60.0, 1.0 / (60.0 * 60.0), //
24.0, 1.0 / 24.0, //
24.0 * 60.0, 1.0 / (24.0 * 60.0), //
24.0 * 60.0 * 60.0, 1.0 / (24.0 * 60.0 * 60.0), };

const CONVERSION_CONSTANTS CONVERSION_CONSTANT = { 180.0 / M_PI, M_PI / 180.0, //
		15.0 * M_PI / 180.0, 1.0 / (15.0 * M_PI / 180.0), //
		15.0, 1.0 / 15.0, //
		15.0 * M_PI / 180.0 / 60.0, 1.0 / (15.0 * M_PI / 180.0 / 60.0), //
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
	return (hour * TIME_CONVERSION_CONSTANT.HOUR_TO_SECOND
			+ minute * TIME_CONVERSION_CONSTANT.MINUTE_TO_SECOND + second)
			* CONVERSION_CONSTANT.SECOND_TO_RADIAN;
}

// trigonometric functions

const double epsilon = 1.0e-14;

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
