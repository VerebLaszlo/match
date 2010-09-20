/*
 * new_match.c
 *
 *  Created on: Sep 18, 2010
 *      Author: vereb
 */

#define restrict __restrict__
#include <stdio.h>
#include <lal/LALDatatypes.h>		// LALStatus
#include <lal/GenerateInspiral.h>	// LALGenerateInspiral()
#include <lal/SimulateCoherentGW.h>	// CoherentGW
#include <lal/LIGOMetadataTables.h>		// SimInspiralTable
#include <lal/GeneratePPNInspiral.h>	// PPNParamStruc
#include <lal/LALNoiseModelsInspiral.h>	// RandomInspiralSignalIn
#include <lal/LALSQTPNWaveformInterface.h>	// XLALSQTPNDestroyCoherentGW()
#include "match.h"

#define PREC "% -15.10lg "
#define PARAMETER_NUMBER 14
#define LINE_MAX 200
#define NUM 12	// number of the predefined configurations


double pi = LAL_PI;

typedef struct tagMatchStruct {
	double *signal[2];
	fftw_complex *csignal[2];
	fftw_plan plan[2];
	int length;
} matchStruct;

typedef struct Parameters {
    int count;
    SimInspiralTable *injParams;
    PPNParamStruc *ppnParams;
} Parameters;

void mallocMatchStruct(matchStruct* m, int length) {
	m->length = length;
	m->signal[0] = fftw_malloc(m->length * sizeof(double));
	m->csignal[0] = fftw_malloc(m->length * sizeof(fftw_complex));
	m->plan[0] = fftw_plan_dft_r2c_1d(m->length, m->signal[0], m->csignal[0],
			FFTW_ESTIMATE);
	m->signal[1] = fftw_malloc(m->length * sizeof(double));
	m->csignal[1] = fftw_malloc(m->length * sizeof(fftw_complex));
	m->plan[1] = fftw_plan_dft_r2c_1d(m->length, m->signal[1], m->csignal[1],
			FFTW_ESTIMATE);
}

void freeMatchStruct(matchStruct *m) {
	fftw_free(m->signal[0]);
	fftw_free(m->csignal[0]);
	fftw_destroy_plan(m->plan[0]);
	fftw_free(m->signal[1]);
	fftw_free(m->csignal[1]);
	fftw_destroy_plan(m->plan[1]);
}

int FillSimInspiralTable(SimInspiralTable *injParams, PPNParamStruc *ppnParams, char **argv) {

        memset(injParams, 0, sizeof(SimInspiralTable));
        memset(ppnParams, 0, sizeof(PPNParamStruc));

	injParams->mass1 = atof(argv[1]);
	injParams->mass2 = atof(argv[2]);
	injParams->spin1x = atof(argv[3]);
	injParams->spin1y = atof(argv[4]);
	injParams->spin1z = atof(argv[5]);
	injParams->spin2x = atof(argv[6]);
	injParams->spin2y = atof(argv[7]);
	injParams->spin2z = atof(argv[8]);
	injParams->qmParameter1 = 1.;
	injParams->qmParameter2 = 1.;
	injParams->inclination = atof(argv[9]);
	injParams->f_lower = atof(argv[10]);
	injParams->distance = atof(argv[11]);
	injParams->polarization = 0;

	ppnParams->deltaT = atof(argv[12]);

        return 0;
}

int ReadFromFile(char *filename , Parameters *params) {
    FILE *datafile;
    char *line;
    // argv[0] is unused
    char *argv[PARAMETER_NUMBER+1];
    int i, counter;

    line = (char*) malloc(sizeof(char)*(LINE_MAX+1));
    datafile = fopen(filename,"r");
    params->count = 0;
    while (fgets(line,LINE_MAX,datafile) != NULL) {
        if (line[0]!='#') {
            params->count++;
        }
    }
    fclose(datafile);

    params->injParams = calloc(sizeof(SimInspiralTable),params->count);
    params->ppnParams = calloc(sizeof(PPNParamStruc),params->count);


    datafile = fopen(filename,"r");
    counter = 0;
    for (i=0; i<PARAMETER_NUMBER; i++) {
        // max 20 digits
        argv[i+1] = malloc(sizeof(char)*20);
    }
    while (fgets(line,LINE_MAX,datafile) != NULL) {
        if (line[0]!='#') {
            sscanf(line,"%s%s%s%s%s%s%s%s%s%s%s%s",
                    argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],argv[8],argv[9],argv[10],argv[11],argv[12]);
            FillSimInspiralTable(&(params->injParams[counter]),&(params->ppnParams[counter]),argv);
            counter++;
        }
    }

    for (i=0; i<PARAMETER_NUMBER; i++) {
        free(argv[i+1]);
    }
    free(line);
    fclose(datafile);

    return 0;
}

int ParseCommandLine(int argc, char **argv, Parameters *params) {
    if (argc<2) {
        printf("Not enough parameter.\n");
        return -1;
    }
    if (!strcmp(argv[1],"--file")) {
        if (argc<3) {
            printf("Please select a file:\n%s --file `filename`\n",argv[0]);
            return -1;
        }
        ReadFromFile(argv[2],params);
    } else { // read from command line
        if (argc<PARAMETER_NUMBER+1) {
            printf("Not enough parameter.\n");
            return -1;
        } else {
            params->count = 1;
            params->injParams = malloc(sizeof(SimInspiralTable));
            params->ppnParams = malloc(sizeof(PPNParamStruc));
            return FillSimInspiralTable(params->injParams, params->ppnParams, argv);
        }
    }

    return 0;
}


int Calculate( SimInspiralTable injParams, PPNParamStruc ppnParams) {
    static LALStatus status;
    double incl = injParams.inclination;
    double chi[2] = { sqrt(SQR(injParams.spin1x) + SQR(injParams.spin1y) + SQR(
                    injParams.spin1z)), sqrt(SQR(injParams.spin2x) + SQR(
                    injParams.spin2y) + SQR(injParams.spin2z)) };

    double freq_i, freq_f, freq_step, fr, df, dt;
    double ts, tt, ss, *psd;
    long index_i, index_f;
    CoherentGW waveform[2];
    static RandomInspiralSignalIn rand;
    matchStruct match;
    FILE *file;
    char filename[50];
    char *PNString[] = { "SpinQuadTaylortwoPointFivePNSS",
                    "SpinQuadTaylortwoPointFivePNALL" };

    freq_i = injParams.f_lower;
    dt = ppnParams.deltaT;
    // antenna függvények kiszámolása
    double t, p, s, fp, fc;
    t = p = s = 0.;
    fp = 0.5 * (1 + t * t) * cos(p) * cos(s) - t * sin(p) * sin(s);
    fc = 0.5 * (1 + t * t) * cos(p) * sin(s) + t * sin(p) * cos(s);
    double hp, hc, a1, a2, phi, shift;
    unsigned i, length;
    short index, second, longer;

    double where[NUM][2][3] = { { { injParams.spin1x / chi[0], injParams.spin1y
                    / chi[0], injParams.spin1z / chi[0] }, { injParams.spin2x / chi[1],
                    injParams.spin2y / chi[1], injParams.spin2z / chi[1] } },
                    { { sin(incl), 0., cos(incl) }, { sin(incl), 0., cos(incl) } }, //1  L_N ~  S_1 ~  S_2
                    { { -sin(incl), 0., -cos(incl) }, { -sin(incl), 0., -cos(incl) } }, //2 -L_N ~  S_1 ~  S_2
                    { { sin(incl), 0., cos(incl) }, { -sin(incl), 0., -cos(incl) } }, //3  L_N ~  S_1 ~ -S_2
                    { { -sin(incl), 0., -cos(incl) }, { sin(incl), 0., cos(incl) } }, //4  L_N ~ -S_1 ~  S_2
                    { { cos(incl), 0., -sin(incl) }, { cos(incl), 0., -sin(incl) } }, //5         S_1 ~  S_2
                    { { cos(incl), 0., -sin(incl) }, { -cos(incl), 0., sin(incl) } }, //6         S_1 ~ -S_2
                    { { 0., 1., 0. }, { 0., 1., 0. } }, // pályasíkos párhuzamos, egyező irány, 	VERY GOOD
                    { { 0., 1., 0. }, { 0., -1., 0. } }, // pályasíkos párhuzamos, különböző irány
                    { { cos(incl), 0., sin(incl) }, { 0., 1., 0. } }, // pályasíkos merőleges						GOOD
                    { { sin(incl), 0., cos(incl) }, { cos(incl), 0., -sin(incl) } },
                    { { cos(incl), 0., -sin(incl) }, { sin(incl), 0., cos(incl) } }, };

    for (index = 0; index < NUM; index++) {
            injParams.spin1x = chi[0] * where[index][0][0];
            injParams.spin1y = chi[0] * where[index][0][1];
            injParams.spin1z = chi[0] * where[index][0][2];
            injParams.spin2x = chi[1] * where[index][1][0];
            injParams.spin2y = chi[1] * where[index][1][1];
            injParams.spin2z = chi[1] * where[index][1][2];
            // hullámformák gyártása
            for (second = 0; second < 2; second++) {
                    memset(&status, 0, sizeof(LALStatus));
                    memset(&waveform[second], 0, sizeof(CoherentGW));
                    LALSnprintf(injParams.waveform, LIGOMETA_WAVEFORM_MAX
                                    * sizeof(CHAR), PNString[second]);
                    LALGenerateInspiral(&status, &waveform[second], &injParams,
                                    &ppnParams);
                    if (status.statusCode) {
                            fprintf(stderr,
                                            "LALSQTPNWaveformTest: error generating waveform\n");
                            return status.statusCode;
                    }
            }

            // hullámformák átírása kezelhetőbb formába
            longer = waveform[0].f->data->length > waveform[1].f->data->length ? 0
                            : 1;
            rand.psd.length = length = waveform[longer].f->data->length;
            freq_f = waveform[longer].f->data->data[length - 1];
            mallocMatchStruct(&match, length);
            for (second = 0; second < 2; second++) {
                    //****  PRÓBA  ****//
                    if (!second) {
                            sprintf(filename, "wave%s%d.out", "ALL", index);
                    } else {
                            sprintf(filename, "wave%s%d.out", "SS", index);
                    }
                    file = fopen(filename, "w");
                    //****  PRÓBA  ****//
                    memset(match.signal[second], 0, match.length * sizeof(double));
                    for (i = 0; i < waveform[second].f->data->length; i++) {
                            a1 = waveform[second].a->data->data[2 * i];
                            a2 = waveform[second].a->data->data[2 * i + 1];
                            phi = waveform[second].phi->data->data[i]
                                            - waveform[second].phi->data->data[0];
                            shift = waveform[second].shift->data->data[i];
                            hp = a1 * cos(shift) * cos(phi) - a2 * sin(shift) * sin(phi);
                            hc = a1 * sin(shift) * cos(phi) + a2 * cos(shift) * sin(phi);
                            match.signal[second][i] = fp * hp + fc * hc;
                            //****  PRÓBA  ****//
                            fprintf(file, PREC PREC"\n", i * dt, match.signal[second][i]);
                            fflush(file);
                            //****  PRÓBA  ****//
                    }
                    XLALSQTPNDestroyCoherentGW(&waveform[second]);
                    fftw_execute(match.plan[second]);
            }
            //****  PRÓBA  ****//
            fclose(file);
            file =fopen("psd.out", "w");
            //****  PRÓBA  ****//

            // karakterisztikus psd előállítása
            df = 1. / dt / rand.psd.length;
            rand.psd.data = (REAL8*) XLALMalloc(sizeof(REAL8) * rand.psd.length);
            LALNoiseSpectralDensity(&status, &rand.psd, &LALLIGOIPsd, df);
            psd = fftw_malloc(length * sizeof(double));
            for (i = 0; i < length; i++) {
                    psd[i] = rand.psd.data[i];
                    fprintf(file, PREC PREC"\n", i*dt, psd[i]);fflush(file);
            }
            //****  PRÓBA  ****//
            fclose(file);
            //****  PRÓBA  ****//

            // overlap kiszámolása
            fr = 0.;
            freq_step = 1. / (dt * length);
            index_i = 0;
            while (fr < freq_i) {
                    fr += freq_step;
                    ++index_i;
            }
            index_f = index_i;
            while (fr < freq_f) {
                    fr += freq_step;
                    index_f++;
            }
            ts = scalar_freq(match.csignal[0], match.csignal[1], psd, freq_i,
                            freq_f);
            tt = scalar_freq(match.csignal[0], match.csignal[0], psd, freq_i,
                            freq_f);
            ss = scalar_freq(match.csignal[1], match.csignal[1], psd, freq_i,
                            freq_f);
//		printf("see: "PREC PREC PREC ", overlap: "PREC"\n", ts, tt, ss, ts
//				/ sqrt(tt * ss));
            printf("overlap: "PREC"\n", ts / sqrt(tt * ss));
            freeMatchStruct(&match);
            XLALFree(rand.psd.data);
            fftw_free(psd);
    }
    puts("Done.");
    return 0;
}

//****  PRÓBA  ****//
//****  PRÓBA  ****//

int main(int argc, char *argv[]) {

        Parameters *params = malloc(sizeof(Parameters));
        int i;
        
	// kezdeti adatok beolvasása
        if (ParseCommandLine(argc,argv,params)<0) {
            return -1;
        }

        for (i=0; i<params->count; i++) {
            Calculate(params->injParams[i],params->ppnParams[i]);
        }


	return 0;
}
