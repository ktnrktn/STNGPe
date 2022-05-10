#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "stngp.h"
#include "def_struct.h"
#include "ion_ch.h"
#include "neuronGP.h"
#include "randnum.h"

#define E_l -60.0
//#define E_l -50.0

void initGP(stateGP *x, char fname[])
{
  double g, v;
  double gNaf = 50.0;
  double gNap = 0.1;
  double gKv2 = 0.1;
  double gKv3 = 10.0;
  double gKv4f = 2.0;
  double gKv4s = 1.0;
  double gKCNQ = 0.2;
  double gCaH = 0.3 ;
  double gHCN = 0.1;
  double gSK = 0.4;
  double gleak = 0.06;
  double A = 3000;
  double sdGP = 0.0;

  FILE *fp;
  char channelName[20];

  if( (fp = fopen(fname,"r")) != NULL){
    while(fscanf(fp, "%s %lf", channelName, &g) == 2){
      if( strcmp(channelName,"Naf") == 0) gNaf = g;
      if( strcmp(channelName,"NaP") == 0) gNap = g;
      if( strcmp(channelName,"Kv2") == 0) gKv2 = g;
      if( strcmp(channelName,"Kv3") == 0) gKv3 = g;
      if( strcmp(channelName,"Kv4f") == 0) gKv4f = g;
      if( strcmp(channelName,"Kv4s") == 0) gKv4s = g;
      if( strcmp(channelName,"KCNQ") == 0) gKCNQ = g;
      if( strcmp(channelName,"CaH") == 0) gCaH = g;
      if( strcmp(channelName,"HCN") == 0) gHCN = g;
      if( strcmp(channelName,"SK") == 0) gSK = g;
      if( strcmp(channelName,"leak") == 0) gleak = g;
      if( strcmp(channelName,"A") == 0) A = g;
      if( strcmp(channelName,"SD") == 0) sdGP = g;
    }
    fclose(fp);
  }

  //v = E_l;
  v = E_l - 10.0;
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
  (*x).x[11] = mKv4f_inf(v);
  (*x).x[12] = hKv4f_inf(v);
  (*x).x[13] = mKv4s_inf(v);
  (*x).x[14] = hKv4s_inf(v);
  (*x).x[15] = mKCNQ_inf(v);
  (*x).x[16] = mCaH_inf(v);
  (*x).x[17] = mHCN_inf(v);
  (*x).x[18] = mSK_inf(Ca0);
  (*x).x[19] = Ca0;

  (*x).gNaf = gNaf*(1.0 + sdGP*gasdev());
  (*x).gNap = gNap*(1.0 + sdGP*gasdev());
  (*x).gKv2 = gKv2*(1.0 + sdGP*gasdev());
  (*x).gKv3 = gKv3*(1.0 + sdGP*gasdev());
  (*x).gKv4f = gKv4f*(1.0 + sdGP*gasdev());
  (*x).gKv4s = gKv4s*(1.0 + sdGP*gasdev());
  (*x).gKCNQ = gKCNQ*(1.0 + sdGP*gasdev());
  (*x).gCaH = gCaH*(1.0 + sdGP*gasdev());
  (*x).gHCN = gHCN*(1.0 + sdGP*gasdev());
  (*x).gSK = gSK*(1.0+sdGP*gasdev());
  (*x).gleak = gleak*(1.0+sdGP*gasdev());
  (*x).A = A*(1.0+sdGP*gasdev());
  //  fprintf(stdout,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", (*x).x[0], gNaf, gNap, gKv2, gKv3, gKv4f, gKv4s, gKCNQ, gCaH, gHCN, gSK, gleak, A); fflush(stdout);

}


//
void derivsGP(stateGP *x, stateGP *dxdt, 
	      double *gAMPA, double *gNMDA, double *gGABA,
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
    mKv4f = x[i].x[11];
    hKv4s = x[i].x[12];
    mKv4f = x[i].x[13];
    hKv4s = x[i].x[14];
    mKCNQ = x[i].x[15];
    mCaH = x[i].x[16];
    mHCN = x[i].x[17];
    mSK = x[i].x[18];
    Ca0 = x[i].x[19];

    g_Naf = x[i].g_Naf;
    g_NaP = x[i].g_NaP;
    g_Kv2 = x[i].g_Kv2;
    g_Kv3 = x[i].g_Kv3;
    g_Kv4f = x[i].g_Kv4f;
    g_Kv4s = x[i].g_Kv4s;
    g_KCNQ = x[i].g_KCNQ;
    g_CaH = x[i].g_CaH;
    g_HCN = x[i].g_HCN;
    g_SK = x[i].g_SK;
    */   
    // g_syn = x[i].g_syn;

    // I_Ca = I_T_GP(gT, v, p, q, v_Ca_GP(Ca_i)) + I_L_GP(gL, v, c, d, v_Ca_GP(Ca_i));
   
    I_syn = gAMPA[i]*(v-E_AMPA) 
      + gNMDA[i]*(v-E_NMDA)/(1.0+0.28*Mg*exp(-0.062*v))
      + gGABA[i]*(v-E_GABAA);
    
    dxdt[i].x[0] = ( 100000.0*Iapp[i]/x[i].A
		     -I_Naf(x[i].x[1], x[i].x[2], x[i].x[3], v, x[i].gNaf)
		     -I_Nap(x[i].x[4], x[i].x[5], x[i].x[6], v, x[i].gNap)
		     -I_Kv2(x[i].x[7], x[i].x[8], v, x[i].gKv2)
		     -I_Kv3(x[i].x[9], x[i].x[10], v, x[i].gKv3)
		     -I_Kv4f(x[i].x[11], x[i].x[12], v, x[i].gKv4f)
		     -I_Kv4s(x[i].x[13], x[i].x[14], v, x[i].gKv4s)
		     -I_KCNQ(x[i].x[15], v, x[i].gKCNQ)
		     -I_CaH(x[i].x[16], v, x[i].gCaH)
		     -I_HCN(x[i].x[17], v, x[i].gHCN)
		     -I_SK(x[i].x[18], v, x[i].gSK)
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
    dxdt[i].x[11] = ( mKv4f_inf(v) -x[i].x[11] )/tau_mKv4f(v);
    dxdt[i].x[12] = ( hKv4f_inf(v) -x[i].x[12] )/tau_hKv4f(v);
    dxdt[i].x[13] = ( mKv4s_inf(v) -x[i].x[13] )/tau_mKv4s(v);
    dxdt[i].x[14] = ( hKv4s_inf(v) -x[i].x[14] )/tau_hKv4s(v);
    dxdt[i].x[15] = ( mKCNQ_inf(v) -x[i].x[15] )/tau_mKCNQ(v);
    dxdt[i].x[16] = ( mCaH_inf(v) -x[i].x[16] )/tau_mCaH(v);
    dxdt[i].x[17] = ( mHCN_inf(v) -x[i].x[17] )/tau_mHCN(v);
    dxdt[i].x[18] = ( mSK_inf(x[i].x[19]) -x[i].x[18] )/tau_mSK(x[i].x[19]);
    dxdt[i].x[19] = -alpha*a2v*I_CaH(x[i].x[16], v, x[i].gCaH) -K_Ca*(x[i].x[19]-Ca0);
  }
}


stateGP addGP(stateGP x, stateGP y){
  int i;
  stateGP z;

  z = x;
  
  for(i=0;i<NumVarGP;i++)
    z.x[i] = x.x[i]+y.x[i];
  
  return z;
}


stateGP mulGP(double x, stateGP y){
  int i;
  stateGP z;

  z = y;

  for(i=0;i<NumVarGP;i++)
    z.x[i] = x*y.x[i];
  
  return z;
}


//rkstepによる近似計算
void rkStepGP(stateGP *x, stateGP *xout,
	      double *gAMPA, double *gNMDA, double *gGABA, double *Iapp,
	      double t, int n)
{
  int i;
  //  static stateGP dxdt[NumGP], dxm[NumGP], dxt[NumGP], xt[NumGP];
  stateGP *dxdt, *dxm, *dxt, *xt;
  double dt2, dt6;
  
  dt2 = dt/2.0;
  dt6 = dt/6.0;
  
  dxdt = (stateGP *) calloc((size_t) n, (size_t) sizeof (stateGP));
  dxm = (stateGP *) calloc((size_t) n, (size_t) sizeof (stateGP));
  dxt = (stateGP *) calloc((size_t) n, (size_t) sizeof (stateGP));
  xt = (stateGP *) calloc((size_t) n, (size_t) sizeof (stateGP));
  
  /* time: t */
  derivsGP(x, dxdt, gAMPA, gNMDA, gGABA, /*gGABAB,*/ Iapp, n);
  for(i=0;i<n;i++) xt[i] = addGP(x[i], mulGP(dt2, dxdt[i]));

  /* time: t + dt/2 */
  derivsGP(xt, dxt, gAMPA, gNMDA, gGABA, /*gGABAB,*/ Iapp, n);
  for(i=0;i<n;i++) xt[i] = addGP(x[i], mulGP(dt2, dxt[i]));

  /* time: t + dt/2 */
  derivsGP(xt, dxm, gAMPA, gNMDA, gGABA, /*gGABAB,*/ Iapp, n);
  for(i=0;i<n;i++){
    xt[i] = addGP(x[i], mulGP(dt, dxm[i]));
    dxm[i] = addGP(dxm[i], dxt[i]);
  }

  /* time: t + dt */
  derivsGP(xt, dxt, gAMPA, gNMDA, gGABA, /*gGABAB,*/ Iapp, n);
  for(i=0;i<n;i++)
    xout[i] = addGP(x[i], mulGP(dt6, addGP(dxdt[i], addGP(dxt[i], mulGP(2.0, dxm[i])))));

  free(dxdt);
  free(dxm);
  free(dxt);
  free(xt);
}


int spikeJudgeGP(FILE *fp, const int cntr, stateGP gp[], 
		  stateGP gp_out[], int spkGP[], int myRank, int num)
{
  int i, flg = 0;
  
  for(i=0;i<num;i++){

    if( (gp[i].x[0] < Eth)&&(gp_out[i].x[0] >= Eth) ){
      fprintf(fp,"%.2lf %d\n", (double) cntr*dt, i); fflush(fp);
      spkGP[i] = 1;
      flg = 1;
    }else{
      spkGP[i] = 0;
    }

    gp[i] = gp_out[i];
    //    if( (i == INDEX) && (myRank == 3) )
    //      gp[i].x[0] = VCLAMP;
    
  }
  return flg;
}

