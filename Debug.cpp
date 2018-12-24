//
// Created by advanced on 17-2-1.
//

#include "MOPSO.h"

void printIn(int option){
    switch (option){
        case 0: cout << "in loops" << endl; break;
        case 1: cout << "in check similarity" << endl; break;
        case 2: cout << "in updateRep" << endl; break;
        case 3: cout << "in updatePBest" << endl; break;
        default: cout << "in no function"; break;
    }
}

void printOut(int option){
    switch (option) {
        case 0: cout << "out loops" << endl;  break;
        case 1: cout << "out check similarity" << endl; break;
        case 2: cout << "out updateRep" << endl;  break;
        case 3: cout << "out updatePBest" << endl; break;
        default: cout << "out no function"; break;
    }
}