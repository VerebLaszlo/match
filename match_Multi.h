/**
 * @file match_Multi.h
 * @author László Veréb
 * @date 2010.03.26.
 */

#include "variables.h"
#include "generator.h"
#include "match.h"
#include <lal/LALNoiseModelsInspiral.h>
#include <lal/GenerateInspiral.h>
#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/LALSQTPNWaveformInterface.h>

#define MODE "SpinTaylorthreePointFivePN" ///<a
#define PN "threePointFivePN"///<a

/**
 * d
 */
typedef enum {
	SPECIAL = 1, SPHERIC = 2,
} cover;

/**
 * d
 */
double calc_Periods(double *per1, double *per2, signalStruct *signal);

/**
 * d
 */
void multi_Match(program_Params *params, binary_System *act, long num, char dir[50]);

/**
 * d
 */
void destroySTWave(CoherentGW waveform);

/**
 * d
 */
binary_System* fill_TDK(long *length, binary_System *minta, cover cov);
