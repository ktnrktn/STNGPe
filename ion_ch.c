#include <stdio.h>
#include <math.h>
#include "ion_ch.h"

/*
double I_l(double v, double g)
{
  return g * (v - E_l);
}
*/

double I_Naf(double m, double h, double s, double v, double g) 
{
  return g * m * m * m * h * s * (v - E_Na);
}
double I_Nap(double m, double h, double s, double v, double g)
{
  return g * m * m * m * h * s * (v - E_Na);
}
double I_Kv2(double m, double h, double v, double g) 
{
  return g * m * m * m * m * h * (v - E_K);
}
double I_Kv3(double m, double h, double v, double g) 
{
  return g * m * m * m * m * h * (v - E_K);
}
double I_Kv4f(double m, double h, double v, double g)
{
  return g * m * m * m * m * h * (v - E_K);
}
double I_Kv4s(double m, double h, double v, double g)
{
  return g * m * m * m * m * h * (v - E_K);
}
double I_KCNQ(double m, double v,double g) 
{
  return g * m * m * m * m * (v - E_K);
}
double I_CaH(double m, double v,double g) 
{
  return g * m * (v - E_Ca);
} 
double I_CaT(double m, double h, double v, double g)
{
  return g * m * m * h * (v - E_Ca);
}
double I_HCN(double m, double v,double g) 
{
  return g * m * (v - E_cat);
}
double I_SK(double m, double v,double g) 
{
  return g * m * (v - E_K);
}

//Naf
double mNaf_inf(double v) 
{
  return 1.0/(1.0 + exp((theta_mNaf - v)/k_mNaf));
}
double tau_mNaf(double v) 
{
  return tau0_mNaf;
}
double hNaf_inf(double v) 
{
  return 1.0/(1.0 + exp((theta_hNaf - v)/k_hNaf));
}
double tau_hNaf(double v) 
{
  return tau0_hNaf + (tau1_hNaf - tau0_hNaf)/(exp((phi_hNaf - v)/sigma0_hNaf) + exp((phi_hNaf - v)/sigma1_hNaf));
}
double sNaf_inf(double v) 
{
  return min_sNaf + (1.0 - min_sNaf)/(1.0 + exp((theta_sNaf - v)/k_sNaf));
}
double tau_sNaf(double v) 
{
  return tau0_sNaf + (tau1_sNaf - tau0_sNaf)/(exp((phi_sNaf - v)/sigma0_sNaf) + exp((phi_sNaf - v)/sigma1_sNaf));
}

//Nap
double mNap_inf(double v) 
{
  return 1.0/(1.0 + exp((theta_mNap - v)/k_mNap));
}
double tau_mNap(double v) 
{
  return tau0_mNap + (tau1_mNap - tau0_mNap)/(exp((phi_mNap - v)/sigma0_mNap) + exp((phi_mNap - v)/sigma1_mNap));
}
double hNap_inf(double v) 
{
  return min_hNap + (1.0 - min_hNap)/(1.0 + exp((theta_hNap - v)/k_hNap));
}
double tau_hNap(double v) 
{
  return tau0_hNap + (tau1_hNap - tau0_hNap)/(exp((phi_hNap - v)/sigma0_hNap) + exp((phi_hNap - v)/sigma1_hNap));
}
double sNap_inf(double v) 
{
  return 1.0/(1.0 + exp((theta_sNap - v)/k_sNap));
}
double alpha_sNap(double v)
{
  return (Aa_sNap*v + Ba_sNap)/(1 - exp((v + Ba_sNap/Aa_sNap)/Ka_sNap));
}
double beta_sNap(double v)
{
  return (Ab_sNap*v + Bb_sNap)/(1 - exp((v + Bb_sNap/Ab_sNap)/Kb_sNap));
}
double tau_sNap(double v) 
{
  return 1.0/(alpha_sNap(v)+beta_sNap(v));
}

//Kv2
double mKv2_inf(double v) 
{
  return 1.0/(1.0 + exp((theta_mKv2 - v)/k_mKv2));
}
double tau_mKv2(double v) 
{
  return tau0_mKv2 + (tau1_mKv2 - tau0_mKv2)/(exp((phi_mKv2 - v)/sigma0_mKv2) + exp((phi_mKv2 - v)/sigma1_mKv2));
}
double hKv2_inf(double v) 
{
  return min_hKv2 + (1.0 - min_hKv2)/(1.0 + exp((theta_hKv2 - v)/k_hKv2));
}
double tau_hKv2(double v) 
{
  return tau0_hKv2;
}

//Kv3
double mKv3_inf(double v) 
{
  return 1.0/(1.0 + exp((theta_mKv3 - v)/k_mKv3));
}
double tau_mKv3(double v) 
{
  return tau0_mKv3 + (tau1_mKv3 - tau0_mKv3)/(exp((phi_mKv3 - v)/sigma0_mKv3) + exp((phi_mKv3 - v)/sigma1_mKv3));
}
double hKv3_inf(double v) 
{
  return min_hKv3 + (1.0 - min_hKv3)/(1.0 + exp((theta_hKv3 - v)/k_hKv3));
}
double tau_hKv3(double v) 
{
  return tau0_hKv3 + (tau1_hKv3 - tau0_hKv3)/(exp((phi_hKv3 - v)/sigma0_hKv3) + exp((phi_hKv3 - v)/sigma1_hKv3));
}

//Kv4f
double mKv4f_inf(double v)
{
  return 1.0/(1.0 + exp((theta_mKv4f - v)/k_mKv4f));
}
double tau_mKv4f(double v)
{
  return tau0_mKv4f + (tau1_mKv4f - tau0_mKv4f)/(exp((phi_mKv4f - v)/sigma0_mKv4f) + exp((phi_mKv4f - v)/sigma1_mKv4f));
}
double hKv4f_inf(double v)
{
  return 1.0/(1.0 + exp((theta_hKv4f - v)/k_hKv4f));
}
double tau_hKv4f(double v)
{
  return tau0_hKv4f + (tau1_hKv4f - tau0_hKv4f)/(exp((phi_hKv4f - v)/sigma0_hKv4f) + exp((phi_hKv4f - v)/sigma1_hKv4f));
  
}
//Kv4s
double mKv4s_inf(double v)
{
  return 1.0/(1.0 + exp((theta_mKv4s - v)/k_mKv4s));
}
double tau_mKv4s(double v)
{
  return tau0_mKv4s + (tau1_mKv4s - tau0_mKv4s)/(exp((phi_mKv4s - v)/sigma0_mKv4s) + exp((phi_mKv4s - v)/sigma1_mKv4s));
}
double hKv4s_inf(double v)
{
  return 1.0/(1.0 + exp((theta_hKv4s - v)/k_hKv4s));
}
double tau_hKv4s(double v)
{
  return tau0_hKv4s + (tau1_hKv4s - tau0_hKv4s)/(exp((phi_hKv4s - v)/sigma0_hKv4s) + exp((phi_hKv4s - v)/sigma1_hKv4s));
}

//KCNQ
double mKCNQ_inf(double v)
{
  return 1.0/(1.0 + exp((theta_mKCNQ - v)/k_mKCNQ));
}
double tau_mKCNQ(double v)
{
  return tau0_mKCNQ + (tau1_mKCNQ - tau0_mKCNQ)/(exp((phi_mKCNQ - v)/sigma0_mKCNQ) + exp((phi_mKCNQ - v)/sigma1_mKCNQ));
}

//CaH
double mCaH_inf(double v) 
{
  return 1.0/(1.0 + exp((theta_mCaH - v)/k_mCaH));
}
double tau_mCaH(double v) 
{
  return tau0_mCaH;
}

//CaT
double mCaT_inf(double v)
{
  return 1.0/(1.0 + exp((theta_mCaT - v)/k_mCaT));
}
double tau_mCaT(double v)
{
  return tau0_mCaT + tau1_mCaT/(exp((phi0_mCaT - v)/sigma0_mCaT) + exp((phi1_mCaT - v)/sigma1_mCaT));
}
double hCaT_inf(double v)
{
  return 1.0/(1.0 + exp((theta_hCaT - v)/k_hCaT));
}
double tau_hCaT(double v)
{
  return tau0_hCaT + tau1_hCaT/(exp((phi_hCaT - v)/sigma0_hCaT) + exp((phi_hCaT - v)/sigma1_hCaT));
}

//HCN
double mHCN_inf(double v) 
{
  return 1.0/(1.0 + exp((theta_mHCN - v)/k_mHCN));
}
double tau_mHCN(double v) 
{
  double tau;
  tau = tau0_mHCN + (tau1_mHCN - tau0_mHCN)/(exp((phi_mHCN - v)/sigma0_mHCN) + exp((phi_mHCN - v)/sigma1_mHCN));
  if(tau < 0.01)
    return 0.01;
  else
    return tau;
}

//SK
double mSK_inf(double Ca) 
{
  double Can;
  Can = pow(Ca, HCoeff);
  
  return Can/(Can + pow(ECh, HCoeff));
}
double tau_mSK(double Ca) 
{
  if(Ca < Ca_sat)
    return tau1_mSK - Ca*(tau1_mSK - tau0_mSK)/Ca_sat;
  else
    return tau0_mSK;
}
