/**
 * @file lal_wrapper.c
 *
 * @date Apr 9, 2011
 * @author vereb
 */

#include "lal_wrapper.h"

void createPSD(LALParameters *lalparams, double *psd) {
	assert(lalparams);
	lalparams->randIn.psd.length = lalparams->max_Length;
	double df = 1. / lalparams->ppnParams->deltaT / lalparams->randIn.psd.length;
	lalparams->randIn.psd.data = (REAL8*)LALMalloc(sizeof(REAL8) * lalparams->randIn.psd.length);
	LALNoiseSpectralDensity(&lalparams->status, &lalparams->randIn.psd, &LALLIGOIPsd, df);
	for (unsigned long j = 0; j < lalparams->randIn.psd.length; j++) {
		psd[j] = lalparams->randIn.psd.data[j];
	}
}
