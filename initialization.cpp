#include "MOPSO.h"

void initializeParticles(Particle *particle){
    for (int j=0; j<Population; j++){
        //  initialize Velocity
        initializeParticle(particle[j]);
        
    }
}

void initializeParticle(Particle &particle){
    int lenP = particle.sizeOfPosition;
    int lenO = particle.sizeOfAddO_origin;

    particle.Velocity = new double[lenP];
    for (int i=0; i < lenP; i++)
        particle.Velocity[i] = 0;
    particle.sizeOfVelocity = lenP;

    //  old_position
    particle.old_position = new double[lenP];

    //  Best
    structBest tmpBest;                                             // apply a new space for tmpBest to store data
    tmpBest.Position = new double [lenP];
    tmpBest.addO_origin = new double [lenO];
    tmpBest.Cost = new double [objectiveNumber];
    cpyDoubleArray(tmpBest.Position,particle.Position,lenP);
    cpyDoubleArray(tmpBest.addO_origin,particle.addO_origin, lenO);
    cpyDoubleArray(tmpBest.Cost, particle.Cost, objectiveNumber);

    particle.Best = tmpBest;                    // particle's pBest point to new particle
    particle.BestArray.push_back(tmpBest);      // put into list
}