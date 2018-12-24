//
// Created by advanced on 17-1-28.
//

#include "MOPSO.h"

void checkAllParticleSimilarity(Particle *particle){
    pthread_t tids[tidSize];
    int haveRun=0;
    while (haveRun < Population){
        for (int i=0; i<tidSize && haveRun+i<Population; i++){
            pthread_create(&tids[i], NULL, checkOneParticleSimilarity, (void *) &particle[haveRun+i]);
        }
        for (int i=0; i<tidSize && haveRun+i<Population; i++){
            pthread_join(tids[i], NULL);
        }
        haveRun += tidSize;
    }
}

void *checkOneParticleSimilarity(void *p) {
    Particle *particle = (Particle *) p;
    printPdb(particle);
    runTM_score(particle->index);
    double TM_score = getTM_score(particle->index);
    cout << "TM_score=" << TM_score << endl;
    // key step:
    if ( ! particleSimilarity(TM_score) ){
        becomeInitialParticle(*particle);
    }
}

void runTM_score(int index){
    const char *particle1 = catStrStr(inputAddress, "1.pdb");
    const char *particle2 = catStrIntStr(tempFileAddress, "temp", index, ".pdb");
    const char *scoreAnswer  = catStrIntStr(tempFileAddress, "score", index, ".txt");
    const char *command = catStrStr("cd ", TM_scoreAddress, " && ./TMscore ", particle1, " ",
                                    particle2, " > ", scoreAnswer);

    system(command);
}

double getTM_score(int index){
    const char *scoreAnswer = catStrIntStr(tempFileAddress, "score", index, ".txt");
    ifstream TM_score(scoreAnswer);
    char buffer[256];
    char ans[256];
    int cnt = 0;
    while (TM_score >> buffer){
        if (strcmp(buffer, "TM-score") == 0){
            cnt++;
        }
        if (cnt == 3) break;
    }
    TM_score >> buffer;
    TM_score >> ans;

    return atof(ans);
}


bool particleSimilarity(double TM_scrore){
    if (TM_scrore > TM_scoreThreshold) return 1;
    else return 0;
}

void becomeInitialParticle(Particle &particle) {
    inputParticle(particle, particle.index, particle.seq);
}

/*double getTM_score(int index){
    const char *scoreAnswer = catStrIntStr(tempFileAddress, "score", index, ".txt");
    freopen(scoreAnswer, "r", stdin);
    char str[10000];
    int cnt = 0;
    do{
        scanf("%s", str);
        if (strcmp(str, "TM-score") == 0) cnt++;
    }while (cnt < 3);                   // find the target TM-score
    scanf("%s", str);                   //  str : "="

    char TM_score[1000];
    scanf("%s", TM_score);
    removeFile(scoreAnswer);
    return atof(TM_score);
}
// bug: consider not find "TM-score" until the end of file / display strcmp
*/