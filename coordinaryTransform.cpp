#include "MOPSO.h"

//  convert rotation to coordinary
void convertRotationToCoordinary(Particle &particle){
    double preR[4][4];
    initialize_1(preR);
    int numAtom = particle.numAtom;
    for (int itAtom=3; itAtom < numAtom; itAtom++){
        double R[4][4];
        double theta = particle.Position[itAtom-3] - particle.old_position[itAtom-3];
        getR(itAtom, particle.origin, theta, R);  // get the new R
        double newR[4][4];
        matrixProduct(R, preR, newR);                                             //  newR = R*preR
        getNewCoordinary(itAtom, newR, particle.origin);            //  (x,y,z,1)t = newR * (x',y',z',1)t

        cpyOriginTOAddO_origin(itAtom, particle.origin, particle.addO_origin);
        if (itAtom % 3 == 0 || itAtom == numAtom-1){
            getNewCoordinaryForO(itAtom, newR, particle.addO_origin, numAtom);
        }
        cpyMatrix(preR, newR);                                                          //  copy newR to preR , a bug there take me 2 hours...T_T
    }
}
