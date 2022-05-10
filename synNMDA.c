#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "def_struct.h"
#include "synNMDA.h"
#include "randnum.h"
#include "stngp.h"

/* NMDA - kinetic model by Destexhe et al. */
#define alphaNMDA 0.5
#define betaNMDA 0.007
#define OpenCntrMax COUNT

extern double sd_g;

/*******************************************************/

synNMDA **allocMtrxNMDA(int n1, int n2)
{
  int i;
  synNMDA **syn;

  syn = calloc((size_t) n1, (size_t) sizeof (synNMDA *));
  for(i=0;i<n1;i++){
    syn[i] = calloc((size_t) n2, (size_t) sizeof (synNMDA));
  }
  
  return syn;
}

void initSynNMDA(synNMDA **syn, int numConnections[], double g, int num)
{
  int i, j;

  double gasdev();
  double genrand_real3();

  for(i=0;i<num;i++)
    for(j=0;j<numConnections[i];j++){
      syn[i][j].r = genrand_real3();
      syn[i][j].openCntr = 0;
      syn[i][j].g = g*(1.0 + sd_g*gasdev());
    }
  
}

void updateNMDA(FILE *fp, synNMDA **syn, int numConnections[],
		int **spkTimeCntr, double **trnsmit,
		double g[], int num)
{
  int i, j;
  double decay;

  decay = exp(-dt*betaNMDA);

  for(i=0;i<num;i++){
    for(j=0;j<numConnections[i];j++){
      if(spkTimeCntr[i][j]==0){
	syn[i][j].openCntr = OpenCntrMax;
	syn[i][j].rs = exp(-dt*(alphaNMDA*trnsmit[i][j]+betaNMDA));
	syn[i][j].rinf = alphaNMDA*trnsmit[i][j]/(alphaNMDA*trnsmit[i][j]+betaNMDA);
      }
      
      if(syn[i][j].openCntr > 0){
	syn[i][j].r = syn[i][j].rinf*(1.0 - syn[i][j].rs) 
	  + syn[i][j].rs*syn[i][j].r;
	syn[i][j].openCntr--;
      }else{
	syn[i][j].r *= decay;
      }
    }
  }

  for(i=0;i<num;i++){
    g[i] = 0.0;
    for(j=0;j<numConnections[i];j++){
      g[i] += syn[i][j].r*syn[i][j].g;
    }
  }
}
