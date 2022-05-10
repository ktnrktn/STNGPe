#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "randnum.h"

int numProc;

int **allocMatrix(int n1, int n2)
{
  int **mtrx;
  int i;
  
  mtrx = calloc((size_t) n1, (size_t) sizeof (int *));
  for(i=0;i<n1;i++) mtrx[i] = calloc((size_t) n2, (size_t) sizeof (int));
  
  return mtrx;
}


void freeMatrix(int **mtrx, int n1, int n2)
{
  int i;

  for(i=0;i<n1;i++) free(mtrx[i]);
  free(mtrx);
}


double distance(int i, int j, int n1, int n2)
{
  double x1, x2;
  double d1, d2;

  x1 = ((double) i)/n1;
  x2 = ((double) j)/n2;

  d1 = x1 - x2;
  if(d1 < 0)
    d1 += 1.0;

  d2 = x2 - x1;
  if(d2 < 0)
    d2 += 1.0;

  if(d1 < d2)
    return d1;
  else
    return d2;
}


void connect(int **link, int n1, int n2, int n, int K)
{
  int i,j,k,cntr,tmp, randno;
  int *candidate;
  double r;
  FILE *fp;

  // n1 = # of post, n2 = # of pre, n = # of connections

  for(i=0;i<n1;i++){

    //    fprintf(stdout,"done [%d]: ", i); fflush(stdout);
    for(j=0, cntr=0; j<n2; j++){
      link[i][j] = 0;

      r = distance(i,j,n1,n2);

      if( r <= ((double) K)/n2 ) // within the range
	cntr++;
    }
    //    fprintf(stdout,"1(%d) ", cntr); fflush(stdout);

    candidate = calloc((size_t) cntr, (size_t) sizeof(int));
    for(j=0, cntr=0; j<n2; j++){

      r = distance(i,j,n1,n2);

      if( r <= ((double) K)/n2 ){
	candidate[cntr] = j;
	cntr++;
      }
    }
    //    fprintf(stdout,"2(%d) ", cntr); fflush(stdout);

    for(j=(cntr-1); j>0; j--){
      randno = genrand_int32()%(j+1);
      tmp = candidate[j];
      candidate[j] = candidate[randno];
      candidate[randno] = tmp;
    }
    //    fprintf(stdout,"3(%d) ", cntr); fflush(stdout);
    
    if(cntr < n) n = cntr;
    for(j=0; j<n; j++){
	link[i][candidate[j]] = 1;
    }
    //    fprintf(stdout,"4(%d) \n", cntr); fflush(stdout);

    free(candidate);
    
  }
  //  fprintf(stdout,"done\n"); fflush(stdout);


  /*
  fp = fopen("topo_ee","w");
  for(i=0;i<n1;i++){
    for(j=0;j<n2;j++){
      if(link[i][j] == 1){
	x = i%vNe; y = i/vNe;
	if(x%16==0 && x!=0 && y%16==0 && y!=0){
	  fprintf(fp, "%d %d\n",x, y);
	  fprintf(fp,"%d %d\n\n",j%vNe, j/vNe);
	}
      }
    }
  }
  fclose(fp);
  */
}


void output(int **link, int n1, int n2, char *fname)
{
  int i,j,k;
  char fname2[10];
  FILE *fp;
  
  if(numProc == 1){
    fp=fopen(fname,"w");
    for(i=0;i<n1;i++){
      for(j=0;j<n2;j++){
	if(link[i][j]==1){
	  fprintf(fp,"%d %d\n",i,j);
	  //	  stat_e[i][0]++;
	}
      }
    }
    fclose(fp);
    
  }else{

    for(k=0;k<numProc;k++){
      sprintf(fname2,"%s.%d",fname,k);
      fp=fopen(fname2,"w");
      for(i=n1*k/numProc;i<n1*(k+1)/numProc;i++){
	for(j=0;j<n2;j++){
	  if(link[i][j]==1){
	    fprintf(fp,"%d %d\n",i%(n1/numProc),j);
	    //	    stat_e[i][0]++;
	  }
	}
      }
      fclose(fp);
    }
  }
}


int main(int argc, char **argv)
{
  int i,j,k;
  int Na, Np, Ns, Ncnnct_GS, Ncnnct_SG, Ncnnct_GG;
  int n_as, n_ps, n_ap, n_aa, n_pp, n_pa, n_sp, n_sa;
  int K_as, K_ps, K_ap, K_aa, K_pp, K_pa, K_sp, K_sa;
  double b1, b2, b3;
  int **link_as, **link_ps, **link_ap, **link_aa, **link_pp, **link_pa, **link_sp, **link_sa;
  int range1, range2, range3, range4;
  unsigned long seed;
  char *fileName[20];
  FILE *fp;
  
  if(argc != 10){
    printf("Usage: %s arg1 ... arg3\n",argv[0]);
    printf("arg1 = bias of STN->GPe connections [1]\n");
    printf("arg2 = bias of Pro->GPe connections [1]\n");
    printf("arg3 = bias of Ark->GPe connections [1]\n");
    printf("arg4 = range of STN->GPe connections (nonselective[1], selective[0])\n");
    printf("arg5 = range of Ark->GPe connections (nonselective[1], selective[0])\n");
    printf("arg6 = range of Pro->GPe connections (nonselective[1], selective[0])\n");
    printf("arg7 = range of Pro->STN connections (nonselective[1], selective[0])\n");
    printf("arg8 = number of process or no output file [-1]\n");
    printf("arg9 = random seed\n");
    exit(0);
  }

  fp = fopen("cparams","w");
  fprintf(fp,"%s %s %s %s %s %s %s %s %s %s\n",
	  argv[0], argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8], argv[9]);

  b1 = (double) atof(argv[1]); // STN to GPe
  b2 = (double) atof(argv[2]); // pre Ark to Pro
  b3 = (double) atof(argv[3]); // post Ark to Pro
  range1 = (int) atoi(argv[4]); // broad or narrow
  range2 = (int) atoi(argv[5]); // broad or narrow
  range3 = (int) atoi(argv[6]); // broad or narrow
  range4 = (int) atoi(argv[7]); // broad or narrow
  numProc = (int) atoi(argv[8]);
  seed = (unsigned long) atoi(argv[9]);

  Na = 768; // 25% of GPe
  Np = 2304;
  Ns = 1024;
  Ncnnct_GS = 35;
  Ncnnct_SG = 100;
  Ncnnct_GG = 65;
   
  /*
  n_as = (int) 35*4/(3*b1+1);
  n_ps = (int) 35*4*b1/(3*b1+1);

  n_ap = (int) 65*3/(3*b2+1);
  n_pp = (int) 65*3*b2/(3*b2+1);
  
  n_aa = (int) 65/(3*b3+1);
  n_pa = (int) 65*b3/(3*b3+1);
  */

  /* original
  n_as = (int) 608/(15*b1+1);
  n_ps = (int) 608*b1/(15*b1+1);
  */

  /* 2020/12/16  miniature = 0.3 nS */
  //n_as = (int) 608/(54*b1+5);
  //n_ps = (int) 608*b1/(54*b1+5);

  /* 2021/01/05  miniature = 0.2 nS */
  n_as = (int) 608/(81*b1+7);
  n_ps = (int) 608*b1/(81*b1+7);
  
  n_ap = (int) 870/(9*b2+4);
  n_pp = (int) 870*b2/(9*b2+4);
  
  n_aa = (int) 58/(3*b3+1);
  n_pa = (int) 58*b3/(3*b3+1);

  n_sp = 100; 
  n_sa = 0; 

  if(range1 == 1){
    K_as = Ns/2;
    K_ps = Ns/2;
  }else{
    K_as = (int) ceil(n_as/2);
    K_ps = (int) ceil(n_ps/2);
  }

  if(range2 == 1){
    K_ap = Np/2;
    K_pp = Np/2;
  }else{
    K_ap = (int) ceil(n_ap/2);
    K_pp = (int) ceil(n_pp/2);
  }

  if(range3 == 1){
    K_aa = Na/2;
    K_pa = Na/2;
  }else{
    K_aa = (int) ceil(n_aa/2);
    K_pa = (int) ceil(n_pa/2);
  }

  if(range4 == 1){
    K_sp = Np/2;
    K_sa = 0;
  }else{
    K_sp = n_sp/2;
    K_sa = 0;
  }


  /*    
  if(argc != 4){
    printf("Usage: %s arg1 ... arg3\n",argv[0]);
    printf("arg1 = type of connectivity pattern [0,1,2]\n");
    printf("arg2 = number of process or no output file [-1]\n");
    printf("arg3 = random seed\n");
    exit(0);
  }

  if(type == 0){
    // model A
    n_as = 35; // Ncnnct_GS
    n_ps = 35; // Ncnnct_GS

    n_ap = 45; // 0.75*Ncnnct_GG
    n_aa = 15; // 0.25*Ncnnct_GG

    n_pp = 45; // 0.75*Ncnnct_GG
    n_pa = 15; // 0.25*Ncnnct_GG

    n_sp = 75; // 0.75*Ncnnct_SG
    n_sa = 25; // 0.25*Ncnnct_SG

    K_as = Ns/2;
    K_ps = Ns/2;
    
    K_ap = Np/2;
    K_aa = Na/2;

    K_pp = Np/2;
    K_pa = Na/2;

    K_sp = Np/2;
    K_sa = Na/2;

  }else if(type == 1){
    // model B
    n_as = 140; // 4*Ncnnct_GS
    n_ps = 0;

    n_ap = 45; // 0.75*Ncnnct_GG
    n_aa = 15; // 0.25*Ncnnct_GG

    n_pp = 45; // 0.75*Ncnnct_GG
    n_pa = 15; // 0.25*Ncnnct_GG

    n_sp = 100; // Ncnnct_SG
    n_sa = 0;

    K_as = Ns/2;
    K_ps = Ns/2;
    
    K_ap = Np/2;
    K_aa = Na/2;

    K_pp = Np/2;
    K_pa = Na/2;

    K_sp = Np/2;
    K_sa = Na/2;

  }else{
    // model C
    n_as = 140; // 4*Ncnnct_GS
    n_ps = 0;

    n_ap = 45; // 0.75*Ncnnct_GG
    n_aa = 15; // 0.25*Ncnnct_GG

    n_pp = 45; // 0.75*Ncnnct_GG
    n_pa = 15; // 0.25*Ncnnct_GG

    n_sp = 100; // Ncnnct_SG
    n_sa = 0;

    K_as = n_as;
    K_ps = n_ps;
    
    K_ap = n_ap;
    K_aa = n_aa;

    K_pp = n_pp;
    K_pa = n_pa;

    K_sp = n_sp;
    K_sa = n_sa;
  }
*/

  init_genrand(seed);

  link_as = allocMatrix(Na,Ns);
  connect(link_as, Na, Ns, n_as, K_as);
  if(numProc > 0) output(link_as, Na, Ns, "mtrx_as");
  freeMatrix(link_as, Na, Ns);

  link_ps = allocMatrix(Np,Ns);
  connect(link_ps, Np, Ns, n_ps, K_ps);
  if(numProc > 0) output(link_ps, Np, Ns, "mtrx_ps");
  freeMatrix(link_ps, Np, Ns);


  link_ap = allocMatrix(Na,Np);
  connect(link_ap, Na, Np, n_ap, K_ap);
  if(numProc > 0) output(link_ap, Na, Np, "mtrx_ap");
  freeMatrix(link_ap, Na, Np);

  link_pp = allocMatrix(Np,Np);
  connect(link_pp, Np, Np, n_pp, K_pp);
  if(numProc > 0) output(link_pp, Np, Np, "mtrx_pp");
  freeMatrix(link_pp, Np, Np);

  link_sp = allocMatrix(Ns,Np);
  connect(link_sp, Ns, Np, n_sp, K_sp);
  if(numProc > 0) output(link_sp, Ns, Np, "mtrx_sp");
  freeMatrix(link_sp, Ns, Np);


  link_aa = allocMatrix(Na,Na);
  connect(link_aa, Na, Na, n_aa, K_aa);
  if(numProc > 0) output(link_aa, Na, Na, "mtrx_aa");
  freeMatrix(link_aa, Na, Na);

  link_pa = allocMatrix(Np,Na);
  connect(link_pa, Np, Na, n_pa, K_pa);
  if(numProc > 0) output(link_pa, Np, Na, "mtrx_pa");
  freeMatrix(link_pa, Np, Na);

  link_sa = allocMatrix(Ns,Na);
  connect(link_sa, Ns, Na, n_sa, K_sa);
  if(numProc > 0) output(link_sa, Ns, Na, "mtrx_sa");
  freeMatrix(link_sa, Ns, Na);


  fclose(fp);

}
