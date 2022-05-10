#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "stngp.h"
#include "def_struct.h"
#include "ion_ch.h"
#include "neuronSTN.h"
#include "randnum.h"

#define E_l -65.0
//#define E_l -65.28


void initSTN(stateSTN *x, char fname[])
{
  double v, g;
  double gNaf = 50.0;
  double gNap = 0.1;
  double gKv2 = 0.1;
  double gKv3 = 10.0;
  double gCaH = 0.3 ;
  double gCaT = 0.3 ;
  double gHCN = 0.1;
  double gSK = 0.4;
  double gleak = 0.06;
  double A = 3000;
  double sdSTN = 0.0;

  char channelName[20];
  FILE *fp;
  
  if( (fp = fopen(fname,"r")) != NULL){
    while(fscanf(fp, "%s %lf", channelName, &g) == 2){
      if( strcmp(channelName,"Naf") == 0) gNaf = g;
      if( strcmp(channelName,"NaP") == 0) gNap = g;
      if( strcmp(channelName,"Kv2") == 0) gKv2 = g;
      if( strcmp(channelName,"Kv3") == 0) gKv3 = g;
      if( strcmp(channelName,"CaH") == 0) gCaH = g;
      if( strcmp(channelName,"CaT") == 0) gCaT = g;
      if( strcmp(channelName,"HCN") == 0) gHCN = g;
      if( strcmp(channelName,"SK") == 0) gSK = g;
      if( strcmp(channelName,"leak") == 0) gleak = g;
      if( strcmp(channelName,"A") == 0) A = g;
      if( strcmp(channelName,"SD") == 0) sdSTN = g;
    }
    fclose(fp);
  }

  //v = E_l;
  v = E_l-10.0;
  //(*x).x[0] = v;
  (*x).x[0] = v*(1.0 + 0.1*gasdev());
  (*x).x[1] = mNaf_inf(v);
  (*x).x[2] = hNaf_inf(v);
  (*x).x[3] = sNaf_inf(v);
  (*x).x[4] = mNap_inf(v);
  (*x).x[5] = hNap_inf(v);
  (*x).x[6] = sNap_inf(v);
  (*x).x[7] = mKv2_inf(v);
  (*x).x[8] = hKv2_inf(v);
  (*x).x[9] = mKv3_inf(v);
  (*x).x[10] = hKv3_inf(v);
  (*x).x[11] = mCaH_inf(v);
  (*x).x[12] = mCaT_inf(v);
  (*x).x[13] = hCaT_inf(v);
  (*x).x[14] = mHCN_inf(v);
  (*x).x[15] = mSK_inf(Ca0);
  (*x).x[16] = Ca0;

  (*x).gNaf = gNaf*(1.0 + sdSTN*gasdev());
  (*x).gNap = gNap*(1.0 + sdSTN*gasdev());
  (*x).gKv2 = gKv2*(1.0 + sdSTN*gasdev());
  (*x).gKv3 = gKv3*(1.0 + sdSTN*gasdev());
  (*x).gCaH = gCaH*(1.0 + sdSTN*gasdev());
  (*x).gCaT = gCaT*(1.0 + sdSTN*gasdev());
  (*x).gHCN = gHCN*(1.0 + sdSTN*gasdev());
  (*x).gSK = gSK*(1.0+sdSTN*gasdev());
  (*x).gleak = gleak*(1.0+sdSTN*gasdev());
  (*x).A = A*(1.0+sdSTN*gasdev());
}


//
void derivsSTN(stateSTN *x, stateSTN *dxdt, 
	      double *gAMPA, double *gNMDA, double *gGABA, /*double *gGABAB,*/
	      double *Iapp, int num)
{
  int i;
  double v;
  // double v, mNaf, hNaf, sNaf, mNap, hNap, sNap, mKv2, hKv2, mKv3, hKv3, mKv4f, hKv4f, mKv4s, hKv4s, mKCNQ, mCaH, mHCN, mSK ;
  double  I_syn; //I_Ca
  //double g_Naf, g_Nap, g_Kv2, g_Kv3, g_Kv4f, g_Kv4s, g_KCNQ, g_CaH, g_HCN,g_SK;
 
  for(i=0;i<num;i++){
    v = x[i].x[0];
    /*
    mNaf = x[i].x[1];
    hNaf = x[i].x[2];
    sNaf = x[i].x[3];
    mNap = x[i].x[4];
    hNap = x[i].x[5];
    sNap = x[i].x[6];
    mKv2 = x[i].x[7];
    hKv2 = x[i].x[8];
    mKv3 = x[i].x[9];
    hKv3 = x[i].x[10];
    mCaH = x[i].x[11];
    mCaT = x[i].x[12];
    hCaT = x[i].x[13];
    mHCN = x[i].x[14];
    mSK = x[i].x[15];
    Ca0 = x[i].x[16];

    g_Naf = x[i].g_Naf;
    g_NaP = x[i].g_NaP;
    g_Kv2 = x[i].g_Kv2;
    g_Kv3 = x[i].g_Kv3;
    g_CaH = x[i].g_CaH;
    g_CaT = x[i].g_CaT;
    g_HCN = x[i].g_HCN;
    g_SK = x[i].g_SK;
    */   
    // g_syn = x[i].g_syn;

   
    I_syn = gAMPA[i]*(v-E_AMPA) 
      + gNMDA[i]*(v-E_NMDA)/(1.0+0.28*Mg*exp(-0.062*v))
      + gGABA[i]*(v-E_GABAA);
    
    dxdt[i].x[0] = ( 100000.0*Iapp[i]/x[i].A
		     -I_Naf(x[i].x[1], x[i].x[2], x[i].x[3], v, x[i].gNaf)
		     -I_Nap(x[i].x[4], x[i].x[5], x[i].x[6], v, x[i].gNap)
		     -I_Kv2(x[i].x[7], x[i].x[8], v, x[i].gKv2)
		     -I_Kv3(x[i].x[9], x[i].x[10], v, x[i].gKv3)
		     -I_CaH(x[i].x[11], v, x[i].gCaH)
		     -I_CaT(x[i].x[12], x[i].x[13], v, x[i].gCaT)
		     -I_HCN(x[i].x[14], v, x[i].gHCN)
		     -I_SK(x[i].x[15], v, x[i].gSK)
		     -x[i].gleak*(v - E_l)
		     -100.0*I_syn/x[i].A) /Cm;
    
    dxdt[i].x[1] = ( mNaf_inf(v) -x[i].x[1] )/tau_mNaf(v);
    dxdt[i].x[2] = ( hNaf_inf(v) -x[i].x[2] )/tau_hNaf(v);
    dxdt[i].x[3] = ( sNaf_inf(v) -x[i].x[3] )/tau_sNaf(v);
    dxdt[i].x[4] = ( mNap_inf(v) -x[i].x[4] )/tau_mNap(v);
    dxdt[i].x[5] = ( hNap_inf(v) -x[i].x[5] )/tau_hNap(v);
    dxdt[i].x[6] = ( sNap_inf(v) -x[i].x[6] )/tau_sNap(v);
    dxdt[i].x[7] = ( mKv2_inf(v) -x[i].x[7] )/tau_mKv2(v);
    dxdt[i].x[8] = ( hKv2_inf(v) -x[i].x[8] )/tau_hKv2(v);
    dxdt[i].x[9] = ( mKv3_inf(v) -x[i].x[9] )/tau_mKv3(v);
    dxdt[i].x[10] = ( hKv3_inf(v) -x[i].x[10] )/tau_hKv3(v);
    dxdt[i].x[11] = ( mCaH_inf(v) -x[i].x[11] )/tau_mCaH(v);
    dxdt[i].x[12] = ( mCaT_inf(v) -x[i].x[12] )/tau_mCaT(v);
    dxdt[i].x[13] = ( hCaT_inf(v) -x[i].x[13] )/tau_hCaT(v);
    dxdt[i].x[14] = ( mHCN_inf(v) -x[i].x[14] )/tau_mHCN(v);
    dxdt[i].x[15] = ( mSK_inf(x[i].x[16]) -x[i].x[15] )/tau_mSK(x[i].x[16]);
    dxdt[i].x[16] = -alpha*a2v*(I_CaH(x[i].x[11], v, x[i].gCaH)+I_CaT(x[i].x[12], x[i].x[13], v, x[i].gCaT)) -K_Ca*(x[i].x[16]-Ca0);
  }
}


stateSTN addSTN(stateSTN x, stateSTN y){
  int i;
  stateSTN z;

  z = x;
  
  for(i=0;i<NumVarSTN;i++)
    z.x[i] = x.x[i]+y.x[i];
  
  return z;
}


stateSTN mulSTN(double x, stateSTN y){
  int i;
  stateSTN z;

  z = y;

  for(i=0;i<NumVarSTN;i++)
    z.x[i] = x*y.x[i];
  
  return z;
}


//rkstepによる近似計算
void rkStepSTN(stateSTN *x, stateSTN *xout,
	       double *gAMPA, double *gNMDA, double *gGABA, double *Iapp, 
	       double t, int n)
{
  int i;
  stateSTN *dxdt, *dxm, *dxt, *xt;
  double dt2, dt6;
  
  dt2 = dt/2.0;
  dt6 = dt/6.0;
  
  dxdt = (stateSTN *) calloc((size_t) n, (size_t) sizeof (stateSTN));
  dxm = (stateSTN *) calloc((size_t) n, (size_t) sizeof (stateSTN));
  dxt = (stateSTN *) calloc((size_t) n, (size_t) sizeof (stateSTN));
  xt = (stateSTN *) calloc((size_t) n, (size_t) sizeof (stateSTN));
  
  /* time: t */
  derivsSTN(x, dxdt, gAMPA, gNMDA, gGABA, /*gGABAB,*/ Iapp, n);
  for(i=0;i<n;i++) xt[i] = addSTN(x[i], mulSTN(dt2, dxdt[i]));

  /* time: t + dt/2 */
  derivsSTN(xt, dxt, gAMPA, gNMDA, gGABA, /*gGABAB,*/ Iapp, n);
  for(i=0;i<n;i++) xt[i] = addSTN(x[i], mulSTN(dt2, dxt[i]));

  /* time: t + dt/2 */
  derivsSTN(xt, dxm, gAMPA, gNMDA, gGABA, /*gGABAB,*/ Iapp, n);
  for(i=0;i<n;i++){
    xt[i] = addSTN(x[i], mulSTN(dt, dxm[i]));
    dxm[i] = addSTN(dxm[i], dxt[i]);
  }

  /* time: t + dt */
  derivsSTN(xt, dxt, gAMPA, gNMDA, gGABA, /*gGABAB,*/ Iapp, n);
  for(i=0;i<n;i++)
    xout[i] = addSTN(x[i], mulSTN(dt6, addSTN(dxdt[i], addSTN(dxt[i], mulSTN(2.0, dxm[i])))));

  free(dxdt);
  free(dxm);
  free(dxt);
  free(xt);
}


int spikeJudgeSTN(FILE *fp, const int cntr, stateSTN stn[],
		   stateSTN stn_out[], int spkSTN[], int myRank, int num)
{
  int i, flg = 0;
  
  for(i=0;i<num;i++){
    if( (stn[i].x[0] < Eth)&&(stn_out[i].x[0] >= Eth) ){
      fprintf(fp,"%.2lf %d\n", (double) cntr*dt, i); fflush(fp);
      spkSTN[i] = 1;
      flg = 1;
    }else{
      spkSTN[i] = 0;
    }

    stn[i] = stn_out[i];
    //    if( (i == INDEX) && (myRank == 7) )
    //      stn[i].x[0] = VCLAMP;

  }
  return flg;
}

