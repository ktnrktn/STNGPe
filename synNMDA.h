#define Mg 1.0

synNMDA **allocMtrxNMDA(int, int);

void initSynNMDA(synNMDA **, int [], double, int);

void updateNMDA(FILE *, synNMDA **, int [], int **, double **, double [], int);
