#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "def_struct.h"
#include "synAMPA.h"
#include "randnum.h"
#include "stngp.h"

#define perVesicle 0.008

/* properties of EPSP/IPSP */
/* AMPA - 3 state model by Tsodyks&Markram (1997) */

//#define tauAMPA 2.7
#define tauAMPA 3.0
extern double ave_tau_r_e;
extern double ave_tau_u_e;
extern double ave_U_e;
extern double sd;
extern double sd_g;

/*******************************************************/

synAMPA **allocMtrxAMPA(int n1, int n2)
{
  int i;
  synAMPA **syn;

  syn = calloc((size_t) n1, (size_t) sizeof (synAMPA *));
  for(i=0;i<n1;i++){
    syn[i] = calloc((size_t) n2, (size_t) sizeof (synAMPA));
  }

  return syn;
}

void initSynAMPA(synAMPA **syn, int numConnections[], double g, double scl, int num)
{
  int i, j;
  double tau_e = tauAMPA;
  double tau_r, tau_u, dcy_e, a;

  double gasdev(void);
  double genrand_real3(void);

  dcy_e = exp(-dt/tau_e);

  for(i=0;i<num;i++){
    for(j=0;j<numConnections[i];j++){

      tau_r = ave_tau_r_e*(1.0 + sd*gasdev());
      if(tau_r < (tau_e +0.1)) tau_r = (tau_e +0.1);

      tau_u = ave_tau_u_e*(1.0 + sd*gasdev());
      if(tau_u < (tau_e +0.1)) tau_u = (tau_e +0.1);

      syn[i][j].dcy_r = exp(-dt/tau_r);
      syn[i][j].dcy_u = exp(-dt/tau_u);
      syn[i][j].a = (dcy_e -syn[i][j].dcy_r)*tau_e/(tau_r - tau_e);

      syn[i][j].U = ave_U_e*(1.0 + sd*gasdev())*scl;
      if(syn[i][j].U < 0.0) syn[i][j].U = 0.0;
      if(syn[i][j].U > 1.0) syn[i][j].U = 1.0;

      syn[i][j].e = 0.0;
      syn[i][j].r = genrand_real3();
      syn[i][j].u = syn[i][j].U;

      syn[i][j].g = g*(1.0 + sd_g*gasdev());
    }
  }
}


void changeSynAMPA(synAMPA **syn, int numConnections[], double mul, int num)
{
  int i, j;

  for(i=0;i<num;i++){
    for(j=0;j<numConnections[i];j++){

      syn[i][j].U *= mul;
      if(syn[i][j].U > 1.0) syn[i][j].U = 1.0;

    }
  }
}


void updateAMPA(FILE *fp, synAMPA **syn, int numConnections[], 
		int **spkTimeCntr, double **trnsmit, 
		double g[], int num)
{
  int i, j, k, numVesicles, cntr;
  double decay;

  decay = exp(-dt/tauAMPA);

  for(i=0;i<num;i++){
    for(j=0;j<numConnections[i];j++){
      if( spkTimeCntr[i][j] == 0 ){
	trnsmit[i][j] = syn[i][j].r;
	numVesicles = (int) floor(syn[i][j].r/perVesicle);
	syn[i][j].u += syn[i][j].U*(1.0-syn[i][j].u);
	for(k=0,cntr=0; k<numVesicles; k++){
	  if(genrand_real3() < syn[i][j].u) cntr++;
	}
	syn[i][j].e += cntr*perVesicle;
	syn[i][j].r -= cntr*perVesicle;
      }else{
	syn[i][j].u *= syn[i][j].dcy_u;
	syn[i][j].e *= decay;
	syn[i][j].r = syn[i][j].r*syn[i][j].dcy_r + syn[i][j].a*syn[i][j].e
	  + (1.0-syn[i][j].dcy_r);
      }
    }
  }

  if(num == 288){
    fprintf(fp, "%lf %lf ", syn[0][0].e, syn[0][0].r);
  }
  
  for(i=0;i<num;i++){
    g[i] = 0.0;
    for(j=0;j<numConnections[i];j++){
      g[i] += syn[i][j].e*syn[i][j].g;
    }
  }
}
