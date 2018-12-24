//
// Created by advanced on 17-1-16.
//

#include "MOPSO.h"
#include "inputFilter.h"

void multiplyBetterInput() {
    char *seq;
    seq = inputSeq();                               //input seq
    /*for (int i=0; i < multiplyNumber; i++) {
        Particle particle;
        inputParticle(particle, i % inputSize, seq);        // i % inputSize : copy different input particles averagely
        initializeParticle(particle);
        multiplyParticle(particle, Population - multiplyNumber + i);
    }*/


    int qualifiedNumber = 0;                        // have copied how many particle.txt in input fold
    for (int i = 0; i < inputSize; i++){
        int multiAmount = 0, needAmount = 0;
        if (i == 0)
            needAmount = multiplyNumber / 2;
        else if(i!=inputSize-1)
            needAmount = (multiplyNumber / 2)/(inputSize-1);
        else
            needAmount=multiplyNumber-qualifiedNumber;
        multiAmount = needAmount * 10;
        Particle *particle = new Particle[multiAmount];


        for (int j = 0; j < multiAmount; j++){
            inputParticle(particle[j], i, seq);
            initializeParticle(particle[j]);
            multiplyParticle(particle[j], j);
            const char *fileName = catStrIntStr(answerAddress, "originCopy", j, ".pdb");
            printPdb(&particle[j], fileName);
        }

        int *res = inputFilter(i+1, multiAmount, needAmount);
        for (int i = 0; i < needAmount; i++) {
            printCopiedOrigin(particle[res[i]], Population - multiplyNumber + qualifiedNumber);
            qualifiedNumber++;
        }
        // delete file <originCopy>
        for (int j = 0; j < multiAmount; j++){
            const char *fileName = catStrIntStr(answerAddress, "originCopy", j, ".pdb");
            removeFile(fileName);
        }

        delete[]particle;
    }

}

void multiplyParticle(Particle &particle, int repNum){
    int len = particle.sizeOfPosition;
    cpyDoubleArray(particle.old_position, particle.Position, len);
    //  calculate Velocity, like PSO update
    for (int i=0; i<len; i++){
        double r = rand() *0.7/ double(RAND_MAX);
        particle.Velocity[i] = r;
    }
    checkV(particle.Velocity, particle.sizeOfVelocity);  

    //  calculate the new position
    for (int i=0; i<len; i++)  {
        particle.Position[i] += particle.Velocity[i];
    }

    checkP(particle.Position, particle.Velocity, particle.sizeOfPosition);
    
    convertRotationToCoordinary(particle);

}
