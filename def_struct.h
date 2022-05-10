#define NumVarSTN 17

typedef struct stateSTN {
  double x[NumVarSTN];
  double gNaf;
  double gNap;
  double gKv2;
  double gKv3;
  double gCaH;
  double gCaT;
  double gHCN;
  double gSK;
  double gleak;
  double A;
} stateSTN;


#define NumVarGP 20

typedef struct stateGP {
  double x[NumVarGP];
  double gNaf;
  double gNap;
  double gKv2;
  double gKv3;
  double gKv4f;
  double gKv4s;
  double gKCNQ;
  double gCaH;
  double gHCN;
  double gSK;
  double gleak;
  double A;
} stateGP;

typedef struct synAMPA {
  double e;
  double r;
  double u;
  double U;
  double dcy_r;
  double dcy_u;
  double a;
  double g;
} synAMPA;

typedef struct synNMDA {
  double r;
  double rs;
  double rinf;
  double g;
  int openCntr;
} synNMDA;

typedef struct synGABAa {
  double e;
  double r;
  double u;
  double U;
  double dcy_r;
  double dcy_u;
  double a;
  double g;
} synGABAa;

typedef struct synBg {
  double gmax;
  //  double rate;
  //  double ginf;
  //  double sqrtD;
  double g;
} synBg;
