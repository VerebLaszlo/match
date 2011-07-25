/*
 * @file test.h
 *
 * @date Apr 1, 2011
 * @author vereb
 */

#ifndef TEST_H_
#define TEST_H_

#define STRUCT {\
	int EGY;\
	double PI;\
}

typedef struct COnstants STRUCT COnstants;

typedef struct CONSTANTS CONSTANTS;

void set_Values(void);

int get_Int(void);

double get_Double(void);

COnstants get_Constants(void);

#endif /* TEST_H_ */
