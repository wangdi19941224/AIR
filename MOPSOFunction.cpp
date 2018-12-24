#include "MOPSO.h"

void PSOAdaptionForPhi(Particle &particle, myRep &rep, int it, const int i) {
    double w = getInertiaWeight(it-1, MaxIt);       //  get inertia wight
    int len = particle.sizeOfPosition;
    // int h = rouletteWheel(rep);
    int h = getGBest(particle, rep);
    cpyDoubleArray(particle.old_position, particle.Position, len);
    myRep ::iterator a = rep.begin();
    for (int i=0; i<h; i++) a++;

    //if (it == 18) printPdb(a, 100+i);

    //  calculate the new Velocity
    double tmp = acos(-2.0 * currentIteration_Times / times_for_each_play + 1) / PI;    // double tmp = acos(-2.0 * it / MaxIt + 1) / PI;
    double c1 = c1min + (c1max - c1min) * (1 - tmp);
    double c2 = c2max - (c2max - c2min) * (1 - tmp);
    printf("it : %d  c1 : %lf  c2 : %lf\n", it, c1, c2);
    for (int i=0; i<len; i++){
        double r1 = rand() / double(RAND_MAX) , r2 = rand() / double(RAND_MAX) ;
        if(it==0 && i%3 != 1) particle.Velocity[i] = rand() / double(RAND_MAX);
        particle.Velocity[i] = w * particle.Velocity[i] +
                               c1 * r1* (particle.Best.Position[i] - particle.Position[i]) +
                               c2 * r2 * (a->Position[i] - particle.Position[i]);
    }
    //  check velocity range
    checkV(particle.Velocity, len);

    cout << "it : " << it << "print V:" <<   endl;
    for (int i = 0; i < len; i++){
        printf("V: %lf VelMax: %lf\n", particle.Velocity[i], VelMax[i]);
    }
    printf("\n");

    //  calculate the new position
    for (int i=0; i<len; i++)   particle.Position[i] += particle.Velocity[i];
    //  check position range
    checkP(particle.Position, particle.Velocity, len);

}

void getAllParticleCost(Particle * particle){
    pthread_t tids[tidSize];
    int haveRun=0;
    while (haveRun < Population){
        for (int i=0; i<tidSize && haveRun+i<Population; i++){
            pthread_create(&tids[i], NULL, getAParticleCost, (void *) &particle[haveRun+i]);
        }
        for (int i=0; i<tidSize && haveRun+i<Population; i++){
            pthread_join(tids[i], NULL);
        }
        haveRun += tidSize;
    }
    removeFile(energyOutputAddress);
}

void *getAParticleCost(void *p){
    freopen(energyOutputAddress, "w", stdout);
    Particle *particle = (Particle *) p;
    printPdb(particle);
    runScwrl(particle->index);

    static long rosettaTime = 0, calTime = 0, charmmTime = 0, mybinTime = 0;
    time_t costStartTime = time(NULL), costEndTime_r, costEndTime_cal, costEndTime_cha, costEndTime_my;

    particle->Cost[0] = getCost_rosetta(particle->index);
    costEndTime_r = time(NULL);

    particle->Cost[1] = getCost_calRWplus(particle->index);
    costEndTime_cal = time(NULL);

    particle->Cost[2] = getCost_charmm(particle->index);
    costEndTime_cha = time(NULL);
/*
    particle->Cost[2] = getCost_mybin(particle->index);
    costEndTime_my = time(NULL);
*/
    if (particle->index == 1){
        rosettaTime += costEndTime_r - costStartTime;
        calTime += costEndTime_cal - costEndTime_r;
        charmmTime += costEndTime_cha - costEndTime_cal;
        mybinTime += 0;

        printCostTime(rosettaTime, calTime, charmmTime, mybinTime);
    }

    removeTmp(particle->index);

    freopen(logAddress, "a", stdout);
}

/*void getAllParticleCostForTest(Particle * particle){                   // calculate f(x)
    for (int i=0; i<Population; i++)
    {
        int len = particle[i].sizeOfAddO_origin;
        double ans1 = 0, ans2=0;
        for (int j=0; j<len; j+=3){
            double x= particle[i].addO_origin[j];
            double y= particle[i].addO_origin[j+1];
            double z= particle[i].addO_origin[j+2];
            ans1 += dis1(x,y,z);
            ans2 += dis2(x, y, z);
        }
        particle[i].Cost[0] = ans1;
        particle[i].Cost[1] = ans2;
    }
}
*/
bool isDominated(double *cost1, double *cost2){
    bool flag = 1;                                          //  if flag is true then  cost1 is dominated by cost2
    //  which means  cost1 >= cost2, for all the f(x)
    bool equal = 1;
    for (int i=0; i<objectiveNumber; i++){
        if (cost1[i] < cost2[i]){   flag = 0;  equal = 0;   break;      }
        if ( !doubleEqual(cost1[i], cost2[i]) )  equal = 0;
    }
    if (equal)  flag  =  1;                             //   if cost1 = cost2 ,for all the f(x), then cost1 is not
    return flag;                                          //dominated by cost2, which means 2 is not better than 1
}                                                       // but it will generate many same paticles in the rep
                                                        // so I change it to true
void decideDominated(Particle * particle){
    for (int i=0; i<Population; i++)
        particle[i].dominated = 0;

    for (int i=0; i<Population; i++)
        for (int j=0; j<Population; j++) {
            if (i==j) continue;
            if (isDominated(particle[i].Cost, particle[j].Cost))
                particle[i].dominated = 1;
        }
}


void updatePBest(Particle * particle){
    for (int i=0; i<Population; i++){
        int lenP = particle[i].sizeOfPosition;
        int lenO = particle[i].sizeOfAddO_origin;
        list<structBest>::iterator it = particle[i].BestArray.begin();
        list<structBest>::iterator tmpIt;
        bool update = true;               //  update=-1  :  pre pBest is bet

        while (it != particle[i].BestArray.end()){
            bool erase = false;
            if (isDominated(particle[i].Cost, it->Cost)){   // particle[i].Cost > it->Cost, new particle is impossible to be pBest
                update = false;
                break;
            }
            if (isDominated(it->Cost, particle[i].Cost)){   // new particle is better than current pBestArray[i], so we need to erase this pBest
                tmpIt = it;
                erase = true;
            }
            it++;
            if (erase) particle[i].BestArray.erase(tmpIt);
            //if (erase) cout << "erase now !" << endl;
        }

        if (update){                                   // use new particle as pBest and put it into pBestArray
            structBest tmpBest;                                             // apply a new space for tmpBest to store data
            tmpBest.Position = new double [lenP];
            tmpBest.addO_origin = new double [lenO];
            tmpBest.Cost = new double [objectiveNumber];
            cpyDoubleArray(tmpBest.Position,particle[i].Position,lenP);
            cpyDoubleArray(tmpBest.addO_origin,particle[i].addO_origin, lenO);
            cpyDoubleArray(tmpBest.Cost, particle[i].Cost, objectiveNumber);

            particle[i].Best = tmpBest;                    // particle's pBest point to new particle
            particle[i].BestArray.push_back(tmpBest);      // put into list
        }
    }
}
