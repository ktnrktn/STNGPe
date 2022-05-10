void initGP(stateGP *, char []);

void derivsGP(stateGP *, stateGP *, double *, double *, double *, double *, int);

void rkStepGP(stateGP *, stateGP *, double *, double *, double *, double *, double, int);

int spikeJudgeGP(FILE *, const int, stateGP [], stateGP [], int [], int, int);

