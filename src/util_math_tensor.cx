/*
 * @file util_math_tensor.c
 *
 * @date Mar 4, 2011
 * @author vereb
 */

#include "util_math_tensor.h"

inline double angle_Between_Vectors3(double left[], double right[]) {
	return angle_Between_Vectors(left, right, 3);
}

inline double angle_Between_Vectors(double left[], double right[], short size) {
	return acos(dot_Product(left, right, size) / (length_Of(left, size), length_Of(right, size)));
}

inline double unit_Vector3_Of(double vector[], double unit_Vector[]) {
	return unit_Vector_Of(vector, 3, unit_Vector);
}

inline double unit_Vector_Of(double vector[], short size, double unit_Vector[]) {
	double length = length_Of(vector, size);
	normalize_Vector(vector, size, length, unit_Vector);
	return length;
}

inline void normalize_Vector3(double vector[], double normalizing_Constant, double normalized[]) {
	normalize_Vector(vector, 3, normalizing_Constant, normalized);
}

inline void normalize_Vector(double vector[], short size, double normalizing_Constant,
	double normalized[]) {
	assert(normalizing_Constant);
	for (short i = 0; i < size; i++) {
		normalized[i] = vector[i] / normalizing_Constant;
	}
}

inline double length_Of_Vector3(double vector[]) {
	return length_Of(vector, 3);
}

inline double length_Of(double vector[], short size) {
	return sqrt(dot_Product(vector, vector, size));
}

inline double dot_Product3(double left[], double right[]) {
	return dot_Product(left, right, 3);
}

inline double dot_Product(double left[], double right[], short size) {
	double scalar = 0.0;
	for (short i = 0; i < size; i++) {
		scalar += left[i] * right[i];
	}
	return scalar;
}

inline void cross_Product3(double left[], double right[], double result[]) {
	result[0] = left[1] * right[2] - left[2] * right[1];
	result[1] = left[2] * right[0] - left[0] * right[2];
	result[2] = left[0] * right[1] - left[1] * right[0];
}

inline void add_Vectors3(double left[], double right[], double result[]) {
	add_Vectors(left, right, 3, result);
}

inline void add_Vectors(double left[], double right[], short size, double result[]) {
	for (short i = 0; i < size; i++) {
		result[i] = left[i] + right[i];
	}
}

double triple_Product(double vec1[], double vec2[], double vec3[]) {
	return vec1[0] * vec2[1] * vec3[2] + vec1[1] * vec2[2] * vec3[0] + vec1[2] * vec2[0] * vec3[1]
		- vec1[0] * vec2[2] * vec3[1] - vec1[1] * vec2[0] * vec3[2] - vec1[2] * vec2[1] * vec3[0];
}
