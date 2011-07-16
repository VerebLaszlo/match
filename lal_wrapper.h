/**
 * @file lal_wrapper.h
 *
 * @date Apr 9, 2011
 * @author vereb
 */

#ifndef LAL_WRAPPER_H_
#define LAL_WRAPPER_H_

#include <lal/LALSQTPNWaveformInterface.h>
#include <lal/LALNoiseModelsInspiral.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/GenerateInspiral.h>
#include <lal/LALInspiralBank.h>
#include <lal/LALDatatypes.h>
#include <lal/LALInspiral.h>
#include <lal/LALStdlib.h>
#include <lal/RealFFT.h>
#include "datatypes.h"

/** LAL parameters.
 */
typedef struct LALParameters {
	LALStatus status;///<a
	CoherentGW waveform[2];///<a
	SimInspiralTable injParams[2];///<a
	PPNParamStruc ppnParams[1];///<a
	RandomInspiralSignalIn randIn;///<a
	short shorter;///<a
	long min_Length;///<a
	long max_Length;///<a
} LALParameters;

/**
 * Done
 * @param lalparams
 * @param psd
 */
void createPSD(LALParameters *lalparams, double *psd);

/** Sets the signal from CoherentGW's A1A2
 * @param i
 * @param sig
 * @param lal
 * @param F
 */
void setSignal_From_A1A2(short i, signalStruct *sig, LALParameters *lal, double F[]);

/** Sets the signal from CoherentGW's HPHC
 * @param i
 * @param sig
 * @param lal
 * @param F
 */
void setSignal_From_HPHC(short i, signalStruct *sig, LALParameters *lal, double F[]);

#endif /* LAL_WRAPPER_H_ */
