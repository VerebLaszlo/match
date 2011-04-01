/**
 * @file main_Test.c
 *
 * @date Mar 15, 2011
 * @author vereb
 */

#include <stdlib.h>
#include <stdio.h>
#include "util_math.h"
#include "test.h"

int main(int argc, char *argv[]) {
	/*short size = 200;
	double radians[size];
	double cosin[size];
	double sine[size];
	double tangent[size];*/
	puts("Start.");
	/*for (short i = 0; i < size; i++) {
		radians[i] = (double) i * M_PI_4;
		if (radians[i] > 20.0*M_PI) {
			size = i;
			break;
		}
	}
	for (short i = 0; i < size; i++) {
		cosin[i] = cos_good(radians[i]);
		sine[i] = sin_good(radians[i]);
		tangent[i] = tan_good(radians[i]);
	}
	for (short i = 0; i < size; i++) {
		printf("%30.24lg %30.24lg %30.24lg\n", radians[i] / M_PI*180.0, tangent[i], tan(radians[i]));
	}*/
	set_Values();
	COnstants CC = get_Constants();
	printf("%d %lg\n", get_Int(), get_Double());
	printf("%d %lg\n", CC.EGY, CC.PI);
	puts("Done.");
	return EXIT_SUCCESS;
}
