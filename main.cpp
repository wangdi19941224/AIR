#include "MOPSO.h"
#include "DH_rotation.h"
#include "TM_score_Combination_Order.h"


int main(int argc, char *argv[])
{
    cout << "Air Begin ^o^" << endl;
    setTime();
    preDisposeInputParametersAndFiles(argv);
    applyVariable();
    getCombinationOrder();                  // get order for the combination of TM_align
    updateVelocityCheck(0);                 // initialize VelMax[]
    multiplyBetterInput();
    inputParticles(particle);               // Input Data
    updateRep(particle, rep, 0);
    initDH_rotation();

    // MOPSO Main Loop
    for (int it = 1; it <= MaxIt; it++) {
        printCurrentId(it);
        updateVelocityCheck(it-1);                    // update VelMax if necessary

        char*PdbAdress=catStrIntStr(answerAddress,it,"/");
        char make[50];
        if(it%10==1)
        {sprintf(make,"mkdir %s",PdbAdress);
        system(make);}
        for (int i = 0; i < Population; i++) {
            PSOAdaptionForPhi(particle[i], rep, it - 1, i);            //  apply the PSO formula
            //convertRotationToCoordinary(particle[i]);    //  convert rotation to coordinary
            convertRotationToOrigin(particle[i]);           //  convert rotation to coordinate by DH
            char*pdbAddress;
            pdbAddress=catStrIntStr(PdbAdress,i,".pdb");
            if(it%10==1){
            printPdb(particle,pdbAddress);}
        }

        updateRep(particle, rep, it);
        updatePBest(particle);                       //      update pBest
        printTime(1, it);
    }
    outputAnswer(rep);
    printTime(0, 0);
    freeUpSpace(particle);
    return 0;
}

/*
In linux terminate:
      ./AIR /home/advanced/DATA/data_test/TR228 4 3 5
 /*
[argv]<0>       <1>                           <2>                     <3>                 <4>
     ./AIR /home/advanced/Data/TR829    7                            3000                 50
       <input address>                 <standard particle numbers>   <iteration times>    <total population>

 */
