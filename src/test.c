/*
 * @file test.c
 *
 * @date Apr 1, 2011
 * @author vereb
 */

#include "test.h"

struct CONSTANTSSTRUCT;

CONSTANTS XX;

void set_Values(void) {
	XX.EGY = 1;
	XX.PI = 3.14;
}

int get_Int(void) {
	return XX.EGY;
}

double get_Double(void) {
	return XX.PI;
}

COnstants get_Constants(void) {
	COnstants temp;
	temp.EGY = XX.EGY;
	temp.PI = XX.PI;
	return temp;
}
