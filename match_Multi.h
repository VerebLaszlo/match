#include "variables.h"
#include "generator.h"
#include "match.h"
#define restrict __restrict__
#include <lal/LALNoiseModelsInspiral.h>
#include <lal/GenerateInspiral.h>
#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/LALSQTPNWaveformInterface.h>

#define MODE "SpinTaylorthreePointFivePN"
#define PN "threePointFivePN"

double calc_Periods(double *per1, double *per2, signalStruct *signal);

void multi_Match(program_Params *params, binary_System *act, long num, char dir[50]);

void destroySTWave(CoherentGW waveform);
