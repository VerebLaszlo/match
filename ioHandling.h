/**
 * @file ioHandling.h
 * @author László Veréb
 * @date 2010.01.18.
 */

//	TODO befejezni

#ifndef IOHANDLING_H_
#define IOHANDLING_H_

#include "util.h"
#include "generator.h"

/**
 *		Structure containing the programs parameters.
 */
typedef struct {
	long numOfGoodMatch;	///< how many waveforms want we to find with a minimum length and with a minimum match
	long numOfBestMatch;	///< how many best waveforms to be
	double deviation;		///< the deviation of the generated waveform's match from the own match
	double min_time;		///< the generated waveform's minimal length
	double max_time;		///< the generated waveform's maximal length. Not in use.
	double delta_time;		///< the differences in the allocated length of the waveform's memory. Not in use.
	int precision;			///< with this precision will be writed the data. Not in use.
} program_params;

/**
 *		The function reads the parameters needed to run the program from the
 *	"program.init" file.
 * @return	: a structure containig the parameters
 */
program_params init_program(void);

/**
 *		The function reads the parameters needed to generate parameters for
 *		waveform-generation.
 * @param[out]	mode		: the waveform generating method
 * @param[out]	dt			: the sampling time
 * @param[out]	freq_Min	: the lowest frequency
 * @param[out]	freq_Max	: the upper frequency
 */
void init_generator(char *mode, dpc dt, dpc freq_Min, dpc freq_Max);

/**
 *		The functions write the waveform to the corresponding file.
 * @param[in]	index	: index of the generated waveform for the file-name
 * @param[in]	data	: vector containing the waveform
 * @param[in]	l		: the length of the vector
 * @param[in]	dt		: the sampling time
 */
void write_Waveh(long index, double data[], long l, double dt);
void write_Wavel(long index, double data[], long l, double dt);
void write_Wavehx(long index, double data[], long l, double dt);
void write_Wavelx(long index, double data[], long l, double dt);

/**
 *		The functions write the generated data to the corresponding "parametesrs.data"
 *	file.
 * @param[in] i			: index of the generated waveform
 * @param[in] freq_Min	: the lowest frequency
 * @param[in] freq_Max	: the upper frequency
 * @param[in] match		: the match vector
 * @param[in] dev		: the deviation vector
 * @param[in] good		: the good/bad vector
 * @param[in] length	: length of the match, deviation and good/bad vector
 * @param[in] params	: generated parameters
 * @param[in] dt		: sampling time
 */
void write_Gen_Data(ulong i, double freq_Min, double freq_Max, double match[], double
	dev[], int good[], int length, const binary_system * const params, double dt);
void write_Gen_Datax(ulong i, double freq_Min, double freq_Max, double match[], double
	dev[], int good[], int length, const binary_system * const params, double dt);

/**
 *		The functions write the generated data for plotting to the corresponding "to_Plot.data"
 *	file.
 * @param[in] i			: index of the generated waveform
 * @param[in] match		: the match vector
 * @param[in] dev		: the deviation vector
 * @param[in] good		: the good/bad vector
 * @param[in] length	: length of the match, deviation and good/bad vector
 * @param[in] params	: generated parameters
 */
void write_Data_to_Plot(ulong i, double match[], double dev[], int good[], int length, const binary_system * const params);
void write_Data_to_Plotx(ulong i, double match[], double dev[], int good[], int length, const binary_system * const params);

#endif /* IOHANDLING_H_ */
