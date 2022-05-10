#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "def_struct.h"
#include "randnum.h"
#include "synBg.h"
#include "stngp.h"


/* background synaptic inputs - point conducance mode by Destexhe et al. (2001) */
//#define u_e     0.1
#define tau_e   3.0
#define tauR_e  800
#define NumBg_e 100

//#define u_i     0.14
#define tau_i   8.0
#define tauR_i  100
#define tauU_i  700
#define NumBg_i 100

double decay_e, u_e;
double decay_i, u_i;

extern double U_ctx, U_str;
extern double sd, sd_g;


void initSynBgExc(synBg syn[], int num_bg, double rate, double gmax, double DAmod, int num)
{
  int i;
  double expTauR, e_infty, e0;

  u_e = DAmod*U_ctx;
  decay_e = exp(-dt/tau_e);

  expTauR = exp(-1000.0/rate/tauR_e);
  e_infty = u_e*(1.0 - expTauR)/(1 - (1 - u_e)*expTauR);

  e0 = tau_e*num_bg*rate*0.001*e_infty;
    
  for(i=0;i<num;i++){
    syn[i].gmax = gmax;
    syn[i].g = gmax*e0;
  }
}


void initSynBgInh(synBg syn[], int num_bg, double rate, double gmax, double DAmod, int num)
{
  int i;
  double expTauR, expTauU, e_infty, e0;

  u_i = DAmod*U_str;
  decay_i = exp(-dt/tau_i);
  
  expTauR = exp(-1000.0/rate/tauR_i);
  expTauU = exp(-1000.0/rate/tauU_i);
  e_infty = u_i*(1.0 - expTauR)/((1.0 - u_i)*(1.0 - expTauR)*(1.0 - expTauU) + u_i);

  e0 = tau_i*num_bg*rate*0.001*e_infty;
    
  for(i=0;i<num;i++){
    syn[i].gmax = gmax;
    syn[i].g = gmax*e0;
  }
}


void synBgExc(synBg syn[], int num_bg, double rate, int num)
{
  int i;
  double expTauR, e_infty, e0, sqrtD;

  double gasdev();

  expTauR = exp(-1000.0/rate/tauR_e);
  e_infty = u_e*(1.0 - expTauR)/(1 - (1 - u_e)*expTauR);

  e0 = tau_e*num_bg*rate*0.001*e_infty;
  sqrtD = e_infty*sqrt(num_bg*rate*0.001*dt);

  for(i=0;i<num;i++){
    syn[i].g = syn[i].g*decay_e + syn[i].gmax*e0*(1.0 - decay_e) 
      +syn[i].gmax*sqrtD*gasdev();
    if(syn[i].g < 0) syn[i].g = 0.0;
  }
}


void synBgInh(synBg syn[], int num_bg, double rate, int num)
{
  int i;
  double expTauR, expTauU, e_infty, e0, sqrtD;
  
  double gasdev();

  expTauR = exp(-1000.0/rate/tauR_i);
  expTauU = exp(-1000.0/rate/tauU_i);
  e_infty = u_i*(1.0 - expTauR)/((1.0 - u_i)*(1.0 - expTauR)*(1.0 - expTauU) + u_i);

  e0 = tau_i*num_bg*rate*0.001*e_infty;
  sqrtD = e_infty*sqrt(num_bg*rate*0.001*dt);
  
  for(i=0;i<num;i++){
    syn[i].g = syn[i].g*decay_i + syn[i].gmax*e0*(1.0 - decay_i)
      +syn[i].gmax*sqrtD*gasdev();
    if(syn[i].g < 0) syn[i].g = 0.0;
  }
}
