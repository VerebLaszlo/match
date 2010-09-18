/*
*  Copyright (C) 2007 Duncan Brown, Jolien Creighton, Thomas Cokelaer, Michele Vallisneri
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/GenerateInspiral.h>

NRCSID(LALSTPNWaveformTestC, "$Id: LALSTPNWaveformTest.c,v 1.1 2004/05/05 20:06:23 thomas Exp");

int main(int argc,char *argv[])
 {
    static LALStatus    mystatus;
    
    CoherentGW      thewaveform;
    SimInspiralTable    injParams;
    PPNParamStruc       ppnParams;

    const char        *filename = argv[15];//"wave1.dat";
    FILE        *outputfile;
    INT4        i,length;
    REAL8       dt;
    REAL8       a1, a2, phi, shift, phi0;
	
	char PNString [50];

	if(argc!=16) 
	{
		printf("                         1  2  3   4   5   6   7   8   9    10      11      12       13 14      15        \n");
		printf("Correct parameter order: m1 m2 S1x S1y S1z S2x S2y S2z incl f_lower f_final distance dt PNorder outputfile\n");
		return(1);
	} 


    memset( &mystatus, 0, sizeof(LALStatus) );
    memset( &thewaveform, 0, sizeof(CoherentGW) );
    memset( &injParams, 0, sizeof(SimInspiralTable) );
    memset( &ppnParams, 0, sizeof(PPNParamStruc) );

    /* --- first we fill the SimInspiral structure --- */
    
    injParams.mass1 = atof(argv[1]);//10.;
    injParams.mass2 = atof(argv[2]);//10.;

    injParams.spin1x = atof(argv[3]);//0.1;
    injParams.spin1y = atof(argv[4]);//0.2;
    injParams.spin1z = atof(argv[5]);//0.3;
    
    injParams.spin2x = atof(argv[6]);//0.4;
    injParams.spin2y = atof(argv[7]);//0.5;
    injParams.spin2z = atof(argv[8]);//0.6;
    
    injParams.inclination = atof(argv[9]);//0.8;

    /* MV-20060224: I believe this is not used in the SpinTaylor code! */
    injParams.f_lower = atof(argv[10]);//40;
    injParams.f_final = atof(argv[11]);//500.0;
     
    /* this is given in Mpc */    
    injParams.distance = atof(argv[12]);//100.;

    ppnParams.deltaT = atof(argv[13]);//1.0 / 16384;//4096.0;

    injParams.polarization   = 0;
    
	sprintf(PNString,"SpinTaylor%s",argv[14]);
	LALSnprintf(injParams.waveform,LIGOMETA_WAVEFORM_MAX*sizeof(CHAR),PNString);//"SpinTaylortwoPN");

    /* this should not be zero*/
//    injParams.theta0 = atof(argv[15]);//0.01;
//    injParams.phi0   = atof(argv[16]);//0.5;
 
    //fprintf(stderr, "Lower cut-off frequency used will be %fHz\n", injParams.f_lower);

    /* --- now we can call the injection function --- */
    LALGenerateInspiral( &mystatus, &thewaveform, &injParams, &ppnParams );
    if ( mystatus.statusCode )
    {
      fprintf( stderr, "LALSTPNWaveformTest: error generating waveform\n" );
      exit( 1 );
    }

    /* --- and finally save in a file --- */

    outputfile = fopen(filename,"w");

    length  = thewaveform.phi->data->length;
    dt      = thewaveform.phi->deltaT;

    phi0    = thewaveform.phi->data->data[0];

    for(i = 0; i < length; i++) {
        a1  = thewaveform.a->data->data[2*i];
        a2  = thewaveform.a->data->data[2*i+1];
        phi     = thewaveform.phi->data->data[i] - phi0;
        shift   = thewaveform.shift->data->data[i];

        fprintf(outputfile,"%e\t%e\t%e\n",
            i*dt,
            a1*cos(shift)*cos(phi) - a2*sin(shift)*sin(phi),
            a1*sin(shift)*cos(phi) + a2*cos(shift)*sin(phi));
    }

    fclose(outputfile);
    //fprintf(stdout,"waveform saved in %s\n",argv[16]);
    return 0;
}
