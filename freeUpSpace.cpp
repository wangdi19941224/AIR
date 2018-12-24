//
// Created by advanced on 17-2-7.
//

#include "MOPSO.h"

void freeUpSpace(Particle *& particle){
    delete [] particle;
}

void freeUpSpace(double **f, int rowsize){
    for (int i = 0; i < rowsize; i++)
        delete[] f[i];
    delete [] f;
}