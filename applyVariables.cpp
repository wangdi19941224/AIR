//
// Created by advanced on 17-2-16.
//

#include "MOPSO.h"

int *par_firNumber;
int *par_secNumber;
Particle *particle;
ANGLE *angle;
myRep rep;
int *sortAns;
REP_COUNT *rep_count;
double *VelMax;

void applyVariable(){
    particle = new Particle[Population];                       // Population = inputFilesNumber + multiplyFilesNumber
    angle = new ANGLE[nRep * 2  ];
    sortAns = new int[nRep * 2];
    rep_count = new REP_COUNT[nRep * 2];
    angle[0].id = -1;

    VelMax = new double [3 * numAA-3];
}

void releaseSpace(double **p, const int &n){
    for (int i = 0; i < n; i++){
        delete [] p[i];
    }
    delete [] p;
}

void releaseSpace(double *p){
    delete [] p;
}