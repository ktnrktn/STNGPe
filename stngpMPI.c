#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
//#include "mpio.h"
#include "def_struct.h"
#include "neuronSTN.h"
#include "neuronGP.h"
#include "synAMPA.h"
#include "synNMDA.h"
#include "synGABAa.h"
#include "synBg.h"
#include "stngp.h"
#include "randnum.h"

#define PI 3.14159265358979323846

#define SizeOfUnsigned 32 /* bit */
#define VCLAMP (-60.0)

/* time step */
//#define dt DELTA
#define DELAY_COUNT_S2G  120 //3.0msec
#define DELAY_COUNT_G2S  120 //3.0msec
#define DELAY_COUNT_G2G   40 //1.0msec

//#define TRANSIENT 500
#define TRANSIENT 8000
#define INPUTDURATION 100

/* size of network */
#define NumProc 8
#define NumSTN 128
#define NumGlobalSTN NumSTN*NumProc
#define NumGPP 288
#define NumGlobalGPP NumGPP*NumProc
#define NumGPA 96
#define NumGlobalGPA NumGPA*NumProc

/* synaptic properties */
double ave_tau_r_e = 200.0;
double ave_tau_r_i = 17300.0;
double ave_tau_u_e = 100.0;
double ave_tau_u_i = dt;
double ave_U_e = 0.1;
double ave_U_i = 0.1;
double sd = 0.1;
double sd_g = 0.1;

/* background synaptic inputs */
double U_ctx = 0.1;
double U_str = 0.14;

double rateMod_ctx1, rateMod_str1;
double rateAmp_ctx, rateAmp_str;
double rateFreq_ctx, rateFreq_str;

/* DA modulation */
/* synaptic properties */
double sclD1=1.0;
double sclD2=1.0;

/* applied current */
double IappSTN[NumSTN] = {0.0};
double IappGPP[NumGPP] = {0.0};
double IappGPA[NumGPA] = {0.0};
double IappSTN0, IappGP0;
double intensity = 0.0;

/* seed of random number */
unsigned long seed, seedNeuron, seedSynapse;

/*******************************************************/

int **allocMtrxINT(int n1, int n2)
{
  int i;
  int **p;

  p = (int **) calloc((size_t) n1, (size_t) sizeof(int *));
  for(i=0;i<n1;i++) p[i] = (int *) calloc((size_t) n2, (size_t) sizeof(int));

  return p;
}


double **allocMtrxDOUBLE(int n1, int n2)
{
  int i;
  double **p;

  p = (double **) calloc((size_t) n1, (size_t) sizeof(double *));
  for(i=0;i<n1;i++) p[i] = (double *) calloc((size_t) n2, (size_t) sizeof(double));

  return p;
}

/*******************************************************/

int connectionInit(int myRank, char fname0[], int **connections, int numConnections[], int **delay)
{
  int i,j;
  char fname1[10];
  FILE *fp;
  int flg = 0;

  sprintf(fname1,"%s.%d",fname0, myRank);
  if( (fp = fopen(fname1,"r")) != NULL ){
    while(fscanf(fp,"%d %d",&i,&j)!=EOF){
      connections[i][numConnections[i]] = j;
      delay[i][numConnections[i]] = 0;
      numConnections[i]++;
    }
    fclose(fp);
  }else{
    fprintf(stderr,"file not found : %s\n", fname0);
    flg = 1;
  }
  //fprintf(stdout,"[%d] %d\n",myRank,flg); fflush(stdout);
  return flg;
}


void initVariables()
{
  int i;

  for(i=0; i<NumSTN;i++)
    IappSTN[i] = 0.0;

  for(i=0;i<NumGPP;i++)
    IappGPP[i] = 0.0;

  for(i=0;i<NumGPA;i++)
    IappGPA[i] = 0.0;
}


void int2bit(int spk[], int num, unsigned spkCompressGlobal[], int numGlobal)
{
  int i, j;
  unsigned b;

  for(j=0;j<(numGlobal/SizeOfUnsigned);j++) spkCompressGlobal[j]=0;

  for(i=(num-1),j=(num/SizeOfUnsigned-1),b=1;i>=0;i--){
    spkCompressGlobal[j] += spk[i]*b;
    if(i%SizeOfUnsigned != 0){
      b *= 2;
    }else{
      j--;
      b = 1;
    }
  }
}


void bit2int(unsigned spkCompressGlobal[], int spkGlobal[], int numGlobal)
{
  int i, j, k;
  unsigned mask, mask0 = 1 << (SizeOfUnsigned-1);
  
  for(i=0;i<numGlobal/SizeOfUnsigned;i++){
    for(j=0,k=i*SizeOfUnsigned,mask=mask0;j<SizeOfUnsigned;j++,k++){
      spkGlobal[k] = ((spkCompressGlobal[i] & mask) >> (SizeOfUnsigned-1-j));
      mask >>= 1;
    }
  }
}


void spkArrivalCheck(int spkGlobal[], int *spkTimeCntr[NumMaxSynapse],
		     int *connections[NumMaxSynapse], int *numConnections,
		     int *delay[NumMaxSynapse], int n)
{
  int i, j;

  for(i=0;i<n;i++){
    for(j=0;j<numConnections[i];j++){
      if(spkGlobal[connections[i][j]] == 1)
	spkTimeCntr[i][j] = delay[i][j];
      else
	spkTimeCntr[i][j]--;
    }
  }
}


void swa(double *rateMod, double rate, double rateAmp)
{
  int i;

  for(i=0; i<(COUNT*1000); i++){
    if(i<(COUNT*500)){
      rateMod[i] = rate*(1.0 + rateAmp*(1 - 2*exp(-0.02*i/COUNT)));
    }else{
      rateMod[i] = rate*(1.0 - rateAmp*(1 - 2*exp(-0.02*(i-500*COUNT)/COUNT)));
    }
  }

  //rateMod_str[j] = rateMod_tmp[(COUNT*1000+j-COUNT*5)%(COUNT*1000)];
}


void beta(double *rateMod, double rate, double rateAmp, double rateFreq, double phaseDiff)
{
  int i;

  for(i=0; i<(COUNT*1000); i++){
    rateMod[i] = rate*(1.0 + rateAmp*sin(2.0*PI*rateFreq*i/(COUNT*1000)-2.0*PI*phaseDiff));
  }
}


double I_AMPA(double gAMPA, double v, double rev)
{
  return gAMPA*(rev - v);
}

double I_NMDA(double gNMDA, double v, double rev)
{
  return gNMDA*(rev - v)/(1.0 + 0.28*Mg*exp(-0.062*v));
}

double I_GABA(double gGABAA, double v, double rev)
{
  return gGABAA*(rev - v);
}

void output(FILE *fpv, int cntr, stateSTN *stn, double *gAMPA2s, double *gNMDA2s, double *gGABA2s, synBg *scAMPA, stateGP *gpp, double *gAMPA2p, double *gNMDA2p, double *gGABA2p, synBg *psGABAa, stateGP *gpa, double *gAMPA2a, double *gNMDA2a, double *gGABA2a, synBg *asGABAa)
{
  double t, v, p, q;
  double Isyn_s, Isyn_p, Isyn_a;
  double iAMPAs, iNMDAs, iGABAs, ibgs;
  double iAMPAp, iNMDAp, iGABAp, ibgp;
  double iAMPAa, iNMDAa, iGABAa, ibga;
  
  t = dt*cntr;

  //iAMPAs = I_AMPA(gAMPA2s[INDEX],stn[INDEX].x[0], 0.0);
  //iNMDAs = I_NMDA(gNMDA2s[INDEX],stn[INDEX].x[0], 0.0);
  iGABAs = I_GABA(gGABA2s[INDEX],stn[INDEX].x[0], -83.0);
  ibgs = I_AMPA(scAMPA[INDEX].g,stn[INDEX].x[0], 0.0);
  //Isyn_s = iAMPAs + iNMDAs + iGABAs;

  iAMPAp = I_AMPA(gAMPA2p[INDEX],gpp[INDEX].x[0], 0.0);
  iNMDAp = I_NMDA(gNMDA2p[INDEX],gpp[INDEX].x[0], 0.0);
  iGABAp = I_GABA(gGABA2p[INDEX],gpp[INDEX].x[0], -83.0);
  ibgp = I_GABA(psGABAa[INDEX].g,gpp[INDEX].x[0], -83.0);
  //Isyn_p = iAMPAp + iNMDAp + iGABAp;

  iAMPAa = I_AMPA(gAMPA2a[INDEX],gpa[INDEX].x[0], 0.0);
  iNMDAa = I_NMDA(gNMDA2a[INDEX],gpa[INDEX].x[0], 0.0);
  iGABAa = I_GABA(gGABA2a[INDEX],gpa[INDEX].x[0], -83.0);
  ibga = I_GABA(asGABAa[INDEX].g,gpa[INDEX].x[0], -83.0);
  //Isyn_a = iAMPAa + iNMDAa + iGABAa;

  fprintf(fpv,"%.2lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", t,
	  stn[INDEX].x[0], ibgs, iGABAs,
	  gpp[INDEX].x[0], ibgp, iAMPAp, iNMDAp, iGABAp,
  	  gpa[INDEX].x[0], ibga, iAMPAa, iNMDAa, iGABAa);
  fflush(fpv);
}

  
int main(int argc, char **argv)
{
  int myRank;
  int numProc;
  int source;
  MPI_Status status;

  stateSTN stn[NumSTN], stn_out[NumSTN];
  stateGP gpp[NumGPP], gpp_out[NumGPP];
  stateGP gpa[NumGPA], gpa_out[NumGPA];

  double gAMPA2s[NumSTN]={0.0}, gNMDA2s[NumSTN]={0.0};
  double gGABAa2s[NumSTN]={0.0};
  double gAMPA2p[NumGPP]={0.0}, gNMDA2p[NumGPP]={0.0};
  double gGABAa2p[NumGPP]={0.0};
  double gAMPA2a[NumGPA]={0.0}, gNMDA2a[NumGPA]={0.0};
  double gGABAa2a[NumGPA]={0.0};
  double gBgExc2s[NumSTN]={0.0};
  double gBgInh2p[NumGPP]={0.0};
  double gBgInh2a[NumGPA]={0.0};

  synAMPA **psAMPA, **asAMPA;
  synNMDA **psNMDA, **asNMDA;
  synGABAa **apGABAa, **ppGABAa, **spGABAa, **aaGABAa, **paGABAa, **saGABAa;
  synBg scAMPA[NumSTN], pcAMPA[NumGPP], acAMPA[NumGPA];
  synBg psGABAa[NumGPP], asGABAa[NumGPA];
  int num_ctx2s, num_ctx2p, num_ctx2a, num_istr, num_dstr;
  double pcProb[NumGPP]={0.0}, acProb[NumGPA]={0.0};

  int numConnections_sa[NumSTN]={0};
  int numConnections_sp[NumSTN]={0};
  int numConnections_pa[NumGPP]={0};
  int numConnections_pp[NumGPP]={0};
  int numConnections_ps[NumGPP]={0};
  int numConnections_aa[NumGPA]={0};
  int numConnections_ap[NumGPA]={0};
  int numConnections_as[NumGPA]={0};

  int **connections_sa, **connections_sp;
  int **connections_pa, **connections_pp, **connections_ps;
  int **connections_aa, **connections_ap, **connections_as;

  int **spkTimeCntr_sa, **spkTimeCntr_sp;
  int **spkTimeCntr_pa, **spkTimeCntr_pp, **spkTimeCntr_ps;
  int **spkTimeCntr_aa, **spkTimeCntr_ap, **spkTimeCntr_as;

  int **delay_sa, **delay_sp;
  int **delay_pa, **delay_pp, **delay_ps;
  int **delay_aa, **delay_ap, **delay_as;

  double **trnsmit_sa, **trnsmit_sp;
  double **trnsmit_pa, **trnsmit_pp, **trnsmit_ps;
  double **trnsmit_aa, **trnsmit_ap, **trnsmit_as;

  int spkGlobalSTN[NumGlobalSTN]={0};
  int spkGlobalGPP[NumGlobalGPP]={0};
  int spkGlobalGPA[NumGlobalGPA]={0};
  int spkSTN[NumSTN]={0};
  int spkGPP[NumGPP]={0};
  int spkGPA[NumGPA]={0};
  unsigned spkCompressGlobalSTN[NumGlobalSTN/SizeOfUnsigned];
  unsigned spkCompressGlobalGPP[NumGlobalGPP/SizeOfUnsigned];
  unsigned spkCompressGlobalGPA[NumGlobalGPA/SizeOfUnsigned];
  unsigned spkCompressSTN[NumGlobalSTN/SizeOfUnsigned];
  unsigned spkCompressGPP[NumGlobalGPP/SizeOfUnsigned];
  unsigned spkCompressGPA[NumGlobalGPA/SizeOfUnsigned];

  int spkSTN_OR, spkGPP_OR, spkGPA_OR;
  int spkGlobalSTN_OR, spkGlobalGPP_OR, spkGlobalGPA_OR;
  int spkGlobalSTN0[NumGlobalSTN]={0};
  int spkGlobalGPP0[NumGlobalGPP]={0};
  int spkGlobalGPA0[NumGlobalGPA]={0};


  double gAMPA_as, gAMPA_ps, gAMPA_sc, gAMPA_pc, gAMPA_ac;
  double gNMDA_as, gNMDA_ps;
  double gGABAa_pp, gGABAa_ap, gGABAa_pa, gGABAa_aa, gGABAa_sp, gGABAa_ps, gGABAa_as;
  double rate_ctx, rate_istr, rate_dstr;
  double rateMod_ctx[1000*COUNT], rateMod_istr[1000*COUNT], rateMod_dstr[1000*COUNT];
  double rateMod_tmp[1000*COUNT];
  double phaseDiff;

  double t, tmax;
  int cntr,cntrMax;
  int i, j, flg, flgOR;

  FILE *fpv,*fps,*fpp,*fpa,*fpsyn,*fp, *fp2;
  char fname[10];

  void init_genrand(unsigned long);

  /* Initialize MPI */
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProc);


  /* in the case of error */
  if(argc != 34){  
    if(myRank == 0){
      printf("[%d/%d] Usage : %s para1 ... para26\n", myRank, numProc, argv[0]);
      printf("para1 = duration\n");

      printf("para2 = g_AMPA from STN to PGP (1.00 nS) \n");
      printf("para3 = g_AMPA from STN to AGP (0.25 nS) \n");
      printf("para4 = g_AMPA from CTX to STN (3.75? nS) \n");
      printf("para5 = g_AMPA from CTX to PGP (3.75? nS) \n");
      printf("para6 = g_AMPA from CTX to AGP (3.75? nS) \n");

      printf("para7 = g_NMDA from STN to PGP (0.5 nS) \n");
      printf("para8 = g_NMDA from STN to AGP (0.17 nS) \n");

      printf("para9 = g_GABAa from PGP to PGP (4.38 nS) \n");
      printf("para10 = g_GABAa from PGP to AGP (6.67 nS) \n");
      printf("para11 = g_GABAa from AGP to PGP (1.33 nS) \n");
      printf("para12 = g_GABAa from AGP to AGP (1.33 nS) \n");

      printf("para13 = g_GABAa from PGP to STN (6.82 nS) \n");
      printf("para14 = g_GABAa from STR to PGP (3.33 nS) \n");
      printf("para15 = g_GABAa from STR to AGP (0.33 nS) \n");

      printf("para16 = base firing rate of bg inputs from CTX (3.0 spikes/s)\n");
      printf("para17 = base firing rate of bg inputs from iSTR (1.0 spikes/s)\n");
      printf("para18 = base firing rate of bg inputs from dSTR (1.0 spikes/s)\n");

      printf("para19 = number of cortical inputs to STN (100)\n");
      printf("para20 = number of cortical inputs to PGPe (100)\n");
      printf("para21 = number of cortical inputs to AGPe (100)\n");
      printf("para22 = number of striatal inputs to PGPe (100)\n");
      printf("para23 = number of striatal inputs to AGPe (100)\n");

      printf("para24 = amplitude of cortical periodical input (x0.0 - x1.0)\n");
      printf("para25 = freqency of cortical periodical input (20.0 Hz)\n");
      printf("para26 = amplitude of striatal periodical input (x0.0 - x1.0)\n");
      printf("para27 = freqency of striatal periodical input (20.0 Hz)\n");
      printf("para28 = phase difference between ctx and str input (0.0)\n");

      printf("para29 = scale for D1R (1.0)\n");
      printf("para30 = scale for D2R (1.0)\n");

      printf("para31 = seed for synapse parameters\n");
      printf("para32 = seed for initial states of neurons\n");
      printf("para33 = seed for simulation\n");
    }
    MPI_Finalize();
    exit(0);
  }
  if(numProc != NumProc){
    printf("[%d/%d] # of processes is wrong\n", myRank, numProc);
    MPI_Finalize();
    exit(0);
  }

  /* set command arguments 0 */
  if(myRank == 0){
    fp = fopen("params","w");
    fprintf(fp,"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n",
	    argv[0],argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],
	    argv[8],argv[9],argv[10],argv[11],argv[12],argv[13],argv[14],
	    argv[15],argv[16],argv[17],argv[18],argv[19],argv[20],argv[21],
	    argv[22],argv[23],argv[24],argv[25],argv[26],argv[27],argv[28],
	    argv[29],argv[30],argv[31],argv[32],argv[33]);
    fclose(fp);
  }
    
  tmax = (double) atof(argv[1]);
  gAMPA_ps = (double) atof(argv[2]);
  gAMPA_as = (double) atof(argv[3]);
  gAMPA_sc = (double) atof(argv[4]);
  gAMPA_pc = (double) atof(argv[5]);
  gAMPA_ac = (double) atof(argv[6]);
  gNMDA_ps = (double) atof(argv[7]);
  gNMDA_as = (double) atof(argv[8]);

  gGABAa_pp = (double) atof(argv[9]);
  gGABAa_ap = (double) atof(argv[10]);
  gGABAa_pa = (double) atof(argv[11]);
  gGABAa_aa = (double) atof(argv[12]);

  gGABAa_sp = (double) atof(argv[13]);
  gGABAa_ps = (double) atof(argv[14]);
  gGABAa_as = (double) atof(argv[15]);
  
  rate_ctx = (double) atof(argv[16]);
  rate_istr = (double) atof(argv[17]);
  rate_dstr = (double) atof(argv[18]);

  num_ctx2s = (int) atoi(argv[19]);
  num_ctx2p = (int) atoi(argv[20]);
  num_ctx2a = (int) atoi(argv[21]);
  num_istr = (int) atoi(argv[22]);
  num_dstr = (int) atoi(argv[23]);
  
  rateAmp_ctx = (double) atof(argv[24]);
  rateFreq_ctx = (double) atof(argv[25]);
  rateAmp_str = (double) atof(argv[26]);
  rateFreq_str = (double) atof(argv[27]);
  phaseDiff = (double) atof(argv[28]);

  sclD1 = (double) atof(argv[29]);
  sclD2 = (double) atof(argv[30]);
  
  seedSynapse = (unsigned long int) atof(argv[31]);
  seedNeuron = (unsigned long int) atoi(argv[32]);
  seed = (unsigned long int) atoi(argv[33]);
  
  cntrMax = (int) tmax/dt;


  /* read files of connection matrix */

  connections_sa = allocMtrxINT(NumSTN, NumMaxSynapse);
  connections_sp = allocMtrxINT(NumSTN, NumMaxSynapse);
  connections_pa = allocMtrxINT(NumGPP, NumMaxSynapse);
  connections_pp = allocMtrxINT(NumGPP, NumMaxSynapse);
  connections_ps = allocMtrxINT(NumGPP, NumMaxSynapse);
  connections_aa = allocMtrxINT(NumGPA, NumMaxSynapse);
  connections_ap = allocMtrxINT(NumGPA, NumMaxSynapse);
  connections_as = allocMtrxINT(NumGPA, NumMaxSynapse);

  delay_sa = allocMtrxINT(NumSTN, NumMaxSynapse);
  delay_sp = allocMtrxINT(NumSTN, NumMaxSynapse);
  delay_pa = allocMtrxINT(NumGPP, NumMaxSynapse);
  delay_pp = allocMtrxINT(NumGPP, NumMaxSynapse);
  delay_ps = allocMtrxINT(NumGPP, NumMaxSynapse);
  delay_aa = allocMtrxINT(NumGPA, NumMaxSynapse);
  delay_ap = allocMtrxINT(NumGPA, NumMaxSynapse);
  delay_as = allocMtrxINT(NumGPA, NumMaxSynapse);
  
  flg = connectionInit(myRank, "mtrx_sa", connections_sa, numConnections_sa, delay_sa)
    || connectionInit(myRank, "mtrx_sp", connections_sp, numConnections_sp, delay_sp)
    || connectionInit(myRank, "mtrx_pa", connections_pa, numConnections_pa, delay_pa)
    || connectionInit(myRank, "mtrx_pp", connections_pp, numConnections_pp, delay_pp)
    || connectionInit(myRank, "mtrx_ps", connections_ps, numConnections_ps, delay_ps)
    || connectionInit(myRank, "mtrx_aa", connections_aa, numConnections_aa, delay_aa)
    || connectionInit(myRank, "mtrx_ap", connections_ap, numConnections_ap, delay_ap)
    || connectionInit(myRank, "mtrx_as", connections_as, numConnections_as, delay_as);

  MPI_Allreduce(&flg,&flgOR,1,MPI_INT,MPI_LOR,MPI_COMM_WORLD);
  if(flgOR != 0){
    MPI_Finalize();
    exit(0);
  }


  /* initialize variables and parameters */
  initVariables();  
  

  /* initialize variables and parameters of synapses */

  spkTimeCntr_sa = allocMtrxINT(NumSTN, NumMaxSynapse);
  spkTimeCntr_sp = allocMtrxINT(NumSTN, NumMaxSynapse);
  spkTimeCntr_pa = allocMtrxINT(NumGPP, NumMaxSynapse);
  spkTimeCntr_pp = allocMtrxINT(NumGPP, NumMaxSynapse);
  spkTimeCntr_ps = allocMtrxINT(NumGPP, NumMaxSynapse);
  spkTimeCntr_aa = allocMtrxINT(NumGPA, NumMaxSynapse);
  spkTimeCntr_ap = allocMtrxINT(NumGPA, NumMaxSynapse);
  spkTimeCntr_as = allocMtrxINT(NumGPA, NumMaxSynapse);

  trnsmit_sa = allocMtrxDOUBLE(NumSTN, NumMaxSynapse);
  trnsmit_sp = allocMtrxDOUBLE(NumSTN, NumMaxSynapse);
  trnsmit_pp = allocMtrxDOUBLE(NumGPP, NumMaxSynapse);
  trnsmit_pa = allocMtrxDOUBLE(NumGPP, NumMaxSynapse);
  trnsmit_ps = allocMtrxDOUBLE(NumGPP, NumMaxSynapse);
  trnsmit_ap = allocMtrxDOUBLE(NumGPA, NumMaxSynapse);
  trnsmit_aa = allocMtrxDOUBLE(NumGPA, NumMaxSynapse);
  trnsmit_as = allocMtrxDOUBLE(NumGPA, NumMaxSynapse);


  seedSynapse *= (myRank+1);
  init_genrand(seedSynapse);

  psAMPA = allocMtrxAMPA(NumGPP,NumMaxSynapse);
  initSynAMPA(psAMPA, numConnections_ps, gAMPA_ps/ave_U_e, sclD2, NumGPP);

  asAMPA = allocMtrxAMPA(NumGPA,NumMaxSynapse);
  initSynAMPA(asAMPA, numConnections_as, gAMPA_as/ave_U_e, sclD2, NumGPA);

  psNMDA = allocMtrxNMDA(NumGPP,NumMaxSynapse);
  initSynNMDA(psNMDA, numConnections_ps, gNMDA_ps, NumGPP);

  asNMDA = allocMtrxNMDA(NumGPA,NumMaxSynapse);
  initSynNMDA(asNMDA, numConnections_as, gNMDA_as, NumGPA);
  
  spGABAa = allocMtrxGABAa(NumSTN,NumMaxSynapse);
  initSynGABAa(spGABAa, numConnections_sp, gGABAa_sp/ave_U_i, sclD2, NumSTN);

  saGABAa = allocMtrxGABAa(NumSTN,NumMaxSynapse);
  initSynGABAa(saGABAa, numConnections_sa, 0.0, sclD2, NumSTN);

  ppGABAa = allocMtrxGABAa(NumGPP,NumMaxSynapse);
  initSynGABAa(ppGABAa, numConnections_pp, gGABAa_pp/ave_U_i, 1.0, NumGPP);

  apGABAa = allocMtrxGABAa(NumGPA,NumMaxSynapse);
  initSynGABAa(apGABAa, numConnections_ap, gGABAa_ap/ave_U_i, 1.0, NumGPA);

  paGABAa = allocMtrxGABAa(NumGPP,NumMaxSynapse);
  initSynGABAa(paGABAa, numConnections_pa, gGABAa_pa/ave_U_i, 1.0, NumGPP);

  aaGABAa = allocMtrxGABAa(NumGPA,NumMaxSynapse);
  initSynGABAa(aaGABAa, numConnections_aa, gGABAa_aa/ave_U_i, 1.0, NumGPA);


  initSynBgExc(scAMPA, num_ctx2s, rate_ctx, gAMPA_sc/U_ctx, sclD2, NumSTN);

  initSynBgExc(pcAMPA, num_ctx2p, rate_ctx, gAMPA_pc/U_ctx, sclD2, NumGPP);

  initSynBgExc(acAMPA, num_ctx2a, rate_ctx, gAMPA_ac/U_ctx, sclD2, NumGPA);

  initSynBgInh(psGABAa, num_istr, rate_istr, gGABAa_ps/U_str, sclD2, NumGPP);

  initSynBgInh(asGABAa, num_dstr, rate_dstr, gGABAa_as/U_str, sclD1, NumGPA);

  /* initialize variables and parameters of neurons */
  seedNeuron *= (myRank+1);
  init_genrand(seedNeuron);

  for(i=0;i<NumSTN;i++) initSTN(&stn[i],"neurontSTN.conf");

  for(i=0;i<NumGPP;i++) initGP(&gpp[i],"neuronGPP.conf");

  for(i=0;i<NumGPA;i++) initGP(&gpa[i],"neuronGPA.conf");

  // swa
  swa(rateMod_ctx, rate_ctx, rateAmp_ctx);
  swa(rateMod_istr, rate_istr, rateAmp_str);
  swa(rateMod_dstr, rate_dstr, rateAmp_str);
  
  
  /* initialize random numbers of Monte-Carlo simulations */
  seed *= (myRank+1);
  init_genrand(seed);


  /* open data files to write */
  sprintf(fname,"v.%d", myRank);
  fpv = fopen(fname,"w");
  sprintf(fname,"sspk.%d",myRank);
  fps = fopen(fname,"w");
  sprintf(fname,"pspk.%d",myRank);
  fpp = fopen(fname,"w");
  sprintf(fname,"aspk.%d",myRank);
  fpa = fopen(fname,"w");
  sprintf(fname,"syn.%d",myRank);
  fp = fopen(fname,"w");


  /* synchronize all the nodes */
  MPI_Barrier(MPI_COMM_WORLD);

  /* start loop */
  for(cntr=0;cntr<cntrMax;cntr++){
    t = dt*cntr;

    for(i=0;i<NumGPP;i++){
      gAMPA2p[i] *= ((double) (cntr-TRANSIENT*COUNT))/1000/COUNT;
      gNMDA2p[i] *= ((double) (cntr-TRANSIENT*COUNT))/1000/COUNT;
      gGABAa2p[i] *= ((double) (cntr-TRANSIENT*COUNT))/1000/COUNT;
    }
    
    for(i=0;i<NumGPA;i++){
      gAMPA2a[i] *= ((double) (cntr-TRANSIENT*COUNT))/1000/COUNT;
      gNMDA2a[i] *= ((double) (cntr-TRANSIENT*COUNT))/1000/COUNT;
      gGABAa2a[i] *= ((double) (cntr-TRANSIENT*COUNT))/1000/COUNT;
      }

    rkStepSTN(stn, stn_out, gAMPA2s, gNMDA2s, gGABAa2s, IappSTN, t, NumSTN);
    rkStepGP(gpp, gpp_out, gAMPA2p, gNMDA2p, gGABAa2p, IappGPP, t, NumGPP);
    rkStepGP(gpa, gpa_out, gAMPA2a, gNMDA2a, gGABAa2a, IappGPA, t, NumGPA);

    spkSTN_OR = spikeJudgeSTN(fps, cntr, stn, stn_out, spkSTN, myRank, NumSTN);
    spkGPP_OR = spikeJudgeGP(fpp, cntr, gpp, gpp_out, spkGPP, myRank, NumGPP);
    spkGPA_OR = spikeJudgeGP(fpa, cntr, gpa, gpa_out, spkGPA, myRank, NumGPA);

    MPI_Allreduce(&spkSTN_OR,&spkGlobalSTN_OR,1,MPI_INT,MPI_LOR,MPI_COMM_WORLD);
    MPI_Allreduce(&spkGPP_OR,&spkGlobalGPP_OR,1,MPI_INT,MPI_LOR,MPI_COMM_WORLD);
    MPI_Allreduce(&spkGPA_OR,&spkGlobalGPA_OR,1,MPI_INT,MPI_LOR,MPI_COMM_WORLD);

    //    if( cntr%5 == 0 && cntr >= (COUNT*10000) )
    if( cntr%8 == 0 )
      output(fpv, cntr, stn, gAMPA2s, gNMDA2s, gGABAa2s, scAMPA,
	     gpp, gAMPA2p, gNMDA2p, gGABAa2p, psGABAa,
	     gpa, gAMPA2a, gNMDA2a, gGABAa2a, asGABAa);

    
    if(spkGlobalSTN_OR == 1){
      int2bit(spkSTN, NumSTN, spkCompressSTN, NumGlobalSTN);
      MPI_Allgather(spkCompressSTN, NumSTN/SizeOfUnsigned, MPI_UNSIGNED,
		    spkCompressGlobalSTN, NumSTN/SizeOfUnsigned, MPI_UNSIGNED,
		    MPI_COMM_WORLD);
      bit2int(spkCompressGlobalSTN, spkGlobalSTN, NumGlobalSTN);
      spkArrivalCheck(spkGlobalSTN, spkTimeCntr_ps, connections_ps, numConnections_ps, delay_ps, NumGPP);
      spkArrivalCheck(spkGlobalSTN, spkTimeCntr_as, connections_as, numConnections_as, delay_as, NumGPA);
    }else{
      spkArrivalCheck(spkGlobalSTN0, spkTimeCntr_ps, connections_ps, numConnections_ps, delay_ps, NumGPP);
      spkArrivalCheck(spkGlobalSTN0, spkTimeCntr_as, connections_as, numConnections_as, delay_as, NumGPA);
    }

    if(spkGlobalGPP_OR == 1){
      int2bit(spkGPP, NumGPP, spkCompressGPP , NumGlobalGPP);
      MPI_Allgather(spkCompressGPP, NumGPP/SizeOfUnsigned, MPI_UNSIGNED,
		    spkCompressGlobalGPP, NumGPP/SizeOfUnsigned, MPI_UNSIGNED,
		    MPI_COMM_WORLD);
      bit2int(spkCompressGlobalGPP, spkGlobalGPP, NumGlobalGPP);
      spkArrivalCheck(spkGlobalGPP, spkTimeCntr_sp, connections_sp, numConnections_sp, delay_sp, NumSTN);
      spkArrivalCheck(spkGlobalGPP, spkTimeCntr_pp, connections_pp, numConnections_pp, delay_pp, NumGPP);
      spkArrivalCheck(spkGlobalGPP, spkTimeCntr_ap, connections_ap, numConnections_ap, delay_ap, NumGPA);
    }else{
      spkArrivalCheck(spkGlobalGPP0, spkTimeCntr_sp, connections_sp, numConnections_sp, delay_sp, NumSTN);
      spkArrivalCheck(spkGlobalGPP0, spkTimeCntr_pp, connections_pp, numConnections_pp, delay_pp, NumGPP);
      spkArrivalCheck(spkGlobalGPP0, spkTimeCntr_ap, connections_ap, numConnections_ap, delay_ap, NumGPA);
    }

    if(spkGlobalGPA_OR == 1){
      int2bit(spkGPA, NumGPA, spkCompressGPA , NumGlobalGPA);
      MPI_Allgather(spkCompressGPA, NumGPA/SizeOfUnsigned, MPI_UNSIGNED,
		    spkCompressGlobalGPA, NumGPA/SizeOfUnsigned, MPI_UNSIGNED,
		    MPI_COMM_WORLD);
      bit2int(spkCompressGlobalGPA, spkGlobalGPA, NumGlobalGPA);
      spkArrivalCheck(spkGlobalGPA, spkTimeCntr_sa, connections_sa, numConnections_sa, delay_sa, NumSTN);
      spkArrivalCheck(spkGlobalGPA, spkTimeCntr_pa, connections_pa, numConnections_pa, delay_pa, NumGPP);
      spkArrivalCheck(spkGlobalGPA, spkTimeCntr_aa, connections_aa, numConnections_aa, delay_aa, NumGPA);
    }else{
      spkArrivalCheck(spkGlobalGPA0, spkTimeCntr_sa, connections_sa, numConnections_sa, delay_sa, NumSTN);
      spkArrivalCheck(spkGlobalGPA0, spkTimeCntr_pa, connections_pa, numConnections_pa, delay_pa, NumGPP);
      spkArrivalCheck(spkGlobalGPA0, spkTimeCntr_aa, connections_aa, numConnections_aa, delay_aa, NumGPA);
    }


    fprintf(fp,"%lf ", t);
    updateAMPA(fp,psAMPA,numConnections_ps,spkTimeCntr_ps,trnsmit_ps,gAMPA2p,NumGPP);
    updateAMPA(fp,asAMPA,numConnections_as,spkTimeCntr_as,trnsmit_as,gAMPA2a,NumGPA);

    updateNMDA(fp,psNMDA,numConnections_ps,spkTimeCntr_ps,trnsmit_ps,gNMDA2p,NumGPP);
    updateNMDA(fp,asNMDA,numConnections_as,spkTimeCntr_as,trnsmit_as,gNMDA2a,NumGPA);

    updateGABAa(fp,spGABAa,numConnections_sp,spkTimeCntr_sp,trnsmit_sp,gGABAa2s,NumSTN,0);
    updateGABAa(fp,ppGABAa,numConnections_pp,spkTimeCntr_pp,trnsmit_pp,gGABAa2p,NumGPP,0);
    updateGABAa(fp,apGABAa,numConnections_ap,spkTimeCntr_ap,trnsmit_ap,gGABAa2a,NumGPA,0);

    updateGABAa(fp,saGABAa,numConnections_sa,spkTimeCntr_sa,trnsmit_sa,gGABAa2s,NumSTN,1);
    updateGABAa(fp,paGABAa,numConnections_pa,spkTimeCntr_pa,trnsmit_pa,gGABAa2p,NumGPP,1);
    updateGABAa(fp,aaGABAa,numConnections_aa,spkTimeCntr_aa,trnsmit_aa,gGABAa2a,NumGPA,1);
    fprintf(fp,"\n");

    synBgExc(scAMPA, num_ctx2s, rateMod_ctx[cntr%(COUNT*1000)], NumSTN);
    synBgExc(pcAMPA, num_ctx2p, rateMod_ctx[cntr%(COUNT*1000)], NumGPP);
    synBgExc(acAMPA, num_ctx2a, rateMod_ctx[cntr%(COUNT*1000)], NumGPA);
    synBgInh(psGABAa, num_istr, rateMod_istr[cntr%(COUNT*1000)], NumGPP);
    synBgInh(asGABAa, num_dstr, rateMod_dstr[cntr%(COUNT*1000)], NumGPA);

    for(i=0;i<NumSTN;i++){
      gAMPA2s[i] = scAMPA[i].g;
    }
    for(i=0;i<NumGPP;i++){
      gAMPA2p[i] += pcAMPA[i].g;
      gGABAa2p[i] += psGABAa[i].g;
    }
    for(i=0;i<NumGPA;i++){
      gAMPA2a[i] += acAMPA[i].g;
      gGABAa2a[i] += asGABAa[i].g;
    }
  }


  fclose(fpv);
  fclose(fps);
  fclose(fpp);
  fclose(fpa);
  fclose(fp);
  
  /* finish MPI */
  MPI_Finalize();

  exit(0);
}
