void initSTN(stateSTN *, char []);

void derivsSTN(stateSTN *, stateSTN *, double *, double *, double *, double *, int);

void rkStepSTN(stateSTN *, stateSTN *, double *, double *, double *, double *, double, int);

int spikeJudgeSTN(FILE *, const int, stateSTN [], stateSTN [], int [], int, int);
