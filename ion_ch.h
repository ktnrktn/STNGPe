#define Cm              1.0 // 0.024 [F/m2]
#define ga 0.2
//#define g_l 0.068

//reversal potential
#define E_Na         50.0
#define E_K         -90.0
#define E_Ca        130.0
#define E_cat       -30.0
//#define E_l         -60.0

//fast Na
#define theta_mNaf  -39.0
#define k_mNaf        5.0
#define tau0_mNaf     0.028
#define tau1_mNaf     0.028

#define theta_hNaf  -48.0
#define k_hNaf       -2.8
#define tau0_hNaf     0.25
#define tau1_hNaf     4.0
#define phi_hNaf    -43.0
#define sigma0_hNaf  10.0
#define sigma1_hNaf  -5.0

#define min_sNaf      0.15
#define theta_sNaf  -40.0
#define k_sNaf       -5.4
#define tau0_sNaf    10.0
#define tau1_sNaf  1000.0
#define phi_sNaf    -40.0
#define sigma0_sNaf  18.3
#define sigma1_sNaf -10.0

//persistent Na
#define theta_mNap  -57.7
#define k_mNap        5.7
#define tau0_mNap     0.03
#define tau1_mNap     0.146
#define phi_mNap    -42.6
#define sigma0_mNap  14.4
#define sigma1_mNap -14.4

#define min_hNap      0.154
#define theta_hNap  -57.0
#define k_hNap       -4.0
#define tau0_hNap    10.0
#define tau1_hNap    17.0
#define phi_hNap    -34.0
#define sigma0_hNap  26.0
#define sigma1_hNap -31.9

#define theta_sNap  -10.0
#define k_sNap       -4.9
#define Aa_sNap      -2.88e-6
#define Ba_sNap      -4.9e-5
#define Ka_sNap       4.63
#define Ab_sNap       6.94e-6
#define Bb_sNap       4.47e-4
#define Kb_sNap      -2.63

//Kv2
#define theta_mKv2  -33.2
#define k_mKv2        9.1
#define tau0_mKv2     0.1
#define tau1_mKv2     3.0
#define phi_mKv2    -33.2
#define sigma0_mKv2  21.7
#define sigma1_mKv2 -13.9

#define min_hKv2      0.2
#define theta_hKv2  -20.0
#define k_hKv2      -10.0
#define tau0_hKv2  3400.0
#define tau1_hKv2  3400.0

//Kv3
#define theta_mKv3  -26.0
#define k_mKv3        7.8
#define tau0_mKv3     0.1
#define tau1_mKv3    14.0
#define phi_mKv3    -26.0
#define sigma0_mKv3  13.0
#define sigma1_mKv3 -12.0

#define min_hKv3      0.6
#define theta_hKv3  -20.0
#define k_hKv3      -10.0
#define tau0_hKv3     7.0
#define tau1_hKv3    33.0
#define phi_hKv3      0.0
#define sigma0_hKv3  10.0
#define sigma1_hKv3 -10.0

//Kv4f
#define theta_mKv4f  -49.0
#define k_mKv4f       12.5
#define tau0_mKv4f     0.25
#define tau1_mKv4f     7.0
#define phi_mKv4f    -49.0
#define sigma0_mKv4f  29.0
#define sigma1_mKv4f -29.0

#define theta_hKv4f  -83.0
#define k_hKv4f      -10.0
#define tau0_hKv4f     7.0
#define tau1_hKv4f    21.0
#define phi_hKv4f    -83.0
#define sigma0_hKv4f  10.0
#define sigma1_hKv4f -10.0

//Kv4s
#define theta_mKv4s  -49.0
#define k_mKv4s       12.5
#define tau0_mKv4s     0.25
#define tau1_mKv4s     7.0
#define phi_mKv4s    -49.0
#define sigma0_mKv4s  29.0
#define sigma1_mKv4s -29.0

#define theta_hKv4s  -83.0
#define k_hKv4s      -10.0
#define tau0_hKv4s    50.0
#define tau1_hKv4s   121.0
#define phi_hKv4s    -83.0
#define sigma0_hKv4s  10.0
#define sigma1_hKv4s -10.0

//KCNQ
#define theta_mKCNQ  -61.0
#define k_mKCNQ       19.5
#define tau0_mKCNQ     6.7
#define tau1_mKCNQ   100.0
#define phi_mKCNQ    -61.0
#define sigma0_mKCNQ  35.0
#define sigma1_mKCNQ -25.0

//CaH
#define theta_mCaH  -20.0
#define k_mCaH        7.0
#define tau0_mCaH     0.2
#define tau1_mCaH     0.2

//CaT
#define theta_mCaT  -56.0
#define k_mCaT        6.7
#define tau0_mCaT     5.0
#define tau1_mCaT     0.33
#define phi0_mCaT   -27.0
#define phi1_mCaT  -102.0
#define sigma0_mCaT -10.0
#define sigma1_mCaT  15.0
#define theta_hCaT  -85.0
#define k_hCaT       -5.8
#define tau0_hCaT     0.0
#define tau1_hCaT   400.0
#define phi_hCaT    -50.0
#define sigma0_hCaT -15.0
#define sigma1_hCaT  16.0

//HCN
#define theta_mHCN  -76.4
#define k_mHCN       -3.3
#define tau0_mHCN     0.0
#define tau1_mHCN  3625.0
#define phi_mHCN    -76.4
#define sigma0_mHCN   6.56
#define sigma1_mHCN  -7.48

//SK
#define ECh           0.35 //[uM]
#define HCoeff        4.6
#define Ca_sat        5.0 //[uM]
#define tau0_mSK      4.0
#define tau1_mSK     76.0

//calcium dynamics
#define Z               2.0
#define F               9.64853415e4//[C mol^{-1}]
#define alpha           1.0 / ( Z * F )
#define a2v          3000.0 //6/d[cm] = 6*10000/d[um]; d=20[um]
//#define a2v          10000.0 //6/d[cm] = 6*10000/d[um]; d=20[um]
#define K_Ca            0.4 // 2.0
#define Ca0             0.01//[uM]


//AMPA NMDA eGABA
#define E_AMPA          0.0
#define E_NMDA          0.0
#define E_GABAA       -83.0
//#define E_syn         0.0
#define Mg              1.0


//Naf
double mNaf_inf(double);
double tau_mNaf(double);

double hNaf_inf(double);
double tau_hNaf(double);

double sNaf_inf(double);
double tau_sNaf(double);


//Nap
double mNap_inf(double);
double tau_mNap(double);

double hNap_inf(double);
double tau_hNap(double);

double sNap_inf(double);
double alpha_sNap(double);
double beta_sNap(double);
double tau_sNap(double);

//Kv2

double mKv2_inf(double);
double tau_mKv2(double);
double hKv2_inf(double);
double tau_hKv2(double);

//Kv3
double mKv3_inf(double);
double tau_mKv3(double);
double hKv3_inf(double);
double tau_hKv3(double);

//Kv4f
double mKv4f_inf(double);
double tau_mKv4f(double);
double hKv4f_inf(double);
double tau_hKv4f(double);

//Kv4s
double mKv4s_inf(double);
double tau_mKv4s(double);
double hKv4s_inf(double);
double tau_hKv4s(double);

//KCNQ
double mKCNQ_inf(double);
double tau_mKCNQ(double);

//CaH
double mCaH_inf(double);
double tau_mCaH(double);

//CaT
double mCaT_inf(double);
double tau_mCaT(double);
double hCaT_inf(double);
double tau_hCaT(double);

//HCN
double mHCN_inf(double);
double tau_mHCN(double);

//SK
double mSK_inf(double);
double tau_mSK(double);


double I_l (double, double);
double I_Naf(double, double, double, double,double);
double I_Nap (double, double, double, double,double);
double I_Kv2 (double, double, double,double);
double I_Kv3 (double, double, double,double);
double I_Kv4f (double, double, double,double);
double I_Kv4s (double, double, double,double);
double I_KCNQ (double, double,double);
double I_CaH (double, double,double);
double I_CaT (double, double, double,double);
double I_HCN (double,double,double);
double I_SK (double, double,double);
