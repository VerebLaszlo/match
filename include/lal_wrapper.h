/**
 * @file lal_wrapper.h
 *
 * @date Apr 9, 2011
 * @author vereb
 */

#ifndef LAL_WRAPPER_H_
#define LAL_WRAPPER_H_

#include "signals.h"
#include "parameters.h"

/**
 * @param parameters
 * @param signal
 * @return
 */
int generateWaveformPair(SystemParameter *parameters, SignalStruct *signal);

#endif /* LAL_WRAPPER_H_ */
