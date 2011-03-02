/**
 * @file match.h
 * @author László Veréb
 * @date 2010.04.08.
 */

#ifndef MATCH_H
#define MATCH_H

#include <fftw3.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/GeneratePPNInspiral.h>
#include "generator.h"

/**	An enum to contains the integer type constatns.
 */
typedef enum {
	H1 = 0,
	H1P = 0,
	H1C = 1,
	H2 = 2,
	H2P = 2,
	H2C = 3,
	HPP = 0,
	HPC = 1,
	HCP = 2,
	HCC = 3,
	NUM_OF_SIGNALS = 4,
} constantValues;

/**	Structure containing the signals.
 */
typedef struct {
	double *signal[NUM_OF_SIGNALS]; ///< the signals in time domain with the following order: \f$h_{1+}\f$, \f$h_{1\times}\f$, \f$h_{2+}\f$, \f$h_{2\times}\f$.
	double *product_Signal[NUM_OF_SIGNALS];///<a
	fftw_complex *csignal[NUM_OF_SIGNALS]; ///< the signals in frequency domain with the following order: \f$\tilde h_{1+}\f$, \f$\tilde h_{1\times}\f$, \f$\tilde h_{2+}\f$, \f$\tilde h_{2\times}\f$.
	fftw_plan plan[NUM_OF_SIGNALS]; ///< FFT plans to calculate \f$\tilde h_{ij}=F\left(h_{ij}\right)\f$
	double *psd;///<d
	long size; ///< the allocated memory for the signals
} signalStruct;

/**
 * Formulae given.
 * \f{gather*}{
 * 	M_{best}=\max_{t_0}\left(\frac{A+B}{2}+\left[\left(\frac{A-B}{2}\right)^2+C^2\right]^{0.5}\right)^{0.5}\;
 * 	M_{minimax}=\max_{t_0}\left(\frac{A+B}{2}-\left[\left(\frac{A-B}{2}\right)^2+C^2\right]^{0.5}\right)^{0.5}\\
 *	A=\inProd{\tilde{e}_{1+}}{\tilde{e}_{2+}}^2+\inProd{\tilde{e}_{1+}}{\tilde{e}_{2\bot}}^2\\
 *	B=\inProd{\tilde{e}_{1\bot}}{\tilde{e}_{2+}}^2+\inProd{\tilde{e}_{1\bot}}{\tilde{e}_{2\bot}}^2\\
 *	C=\inProd{\tilde{e}_{1+}}{\tilde{e}_{2+}}+\inProd{\tilde{e}_{1\bot}}{\tilde{e}_{2+}}+\inProd{\tilde{e}_{1+}}{\tilde{e}_{2\bot}}+\inProd{\tilde{e}_{1\bot}}{\tilde{e}_{2\bot}}
 * \f}
 * @param in
 * @param min_Index
 * @param max_Index
 * @param typ
 * @param best
 * @param minimax
 */
void calc_Matches(signalStruct *in, long min_Index, long max_Index, double *typ, double *best,
		double *minimax);

/**
 * Needs formulae.
 * @param in
 * @param out
 * @param min_Index
 * @param max_Index
 */
void orthonormalise(signalStruct *in, long min_Index, long max_Index, signalStruct *out);

/**
 * Formulae given.
 * \f{gather*}{
 * 	\tilde{e}=\frac{\tilde{h}}{\sqrt{\inProd{\tilde{h}}{\tilde{h}}}}
 * \f}
 * @param in
 * @param out
 * @param min_Index
 * @param max_Index
 */
void normalise(signalStruct *in, long min_Index, long max_Index, signalStruct *out);

/**
 * Formulae given.
 * \f{gather*}{
 * 	\inProd{h_1}{h_2}=4\Re\int_{f_{min}}^{f_{max}}\frac{\tilde{h}_1(f)\tilde{h}_2^*(f)}{S_h(f)}df
 * \f}
 * @param left
 * @param right
 * @param norm
 * @param min_Index
 * @param max_Index
 * @return
 */
double inner_Product(fftw_complex left[], fftw_complex right[], double norm[], long min_Index,
		long max_Index);

/**
 * Needs formulae.
 * \f{gather*}{
 * 	\tilde{e}_\bot=\tilde{e}_\times-\tilde{e}_+\frac{\inProd{\tilde{e}_+}{\tilde{e}_\times}}{\inProd{\tilde{e}_+}{\tilde{e}_+}}
 * \f}
 * @param in
 * @param out
 * @param min_Index
 * @param max_Index
 */
void orthogonise(signalStruct *in, long min_Index, long max_Index, signalStruct *out);

/**
 * Needs formulae.
 * @param left
 * @param right
 * @param norm
 * @param min_Index
 * @param max_Index
 * @param out
 */
void cross_Product(fftw_complex left[], fftw_complex right[], double norm[], long min_Index,
		long max_Index, fftw_complex out[]);

/**
 * Needs formulae.
 * @param in
 * @param min_Index
 * @param max_Index
 * @param typ
 * @param best
 * @param minimax
 */
void calc_Timemaximised_Matches(signalStruct *in, long min_Index, long max_Index, double *typ,
		double *best, double *minimax);

/**	Allocates memory for the signals.
 * @param[in] s		: pointer to the structure containing the signals
 * @param[in] size	: the size of the allocated memory
 * @return : the succes of the memory allocation
 */
void create_Signal_Struct(signalStruct *s, long size);

/**	Deallocates the memory of the signals.
 * @param[in] s	: pointer to the structure containing the signals
 */
void destroy_Signal_Struct(signalStruct *s);

#endif
