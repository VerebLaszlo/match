/*
 * @file util_math.c
 * @date 2010-08-16 10.27.00
 * @author László Veréb
 * @brief Various utilities for math.
 */

#include "util_math.h"

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

inline double random_From_Zero_To(double top) {
	return top * random_From_Zero_To_One();
}

inline double random_Between(double bottom, double top) {
	return (top - bottom) * random_From_Zero_To_One() + bottom;
}

inline double random_From_Zero_To_One(void) {
	return (double) rand() / ((double) RAND_MAX + 1.);
}

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
