#include "MOPSO.h"

//  the address of the input data
const char rootAddress[] = "/home/wangdi/mopso/";
const char *energyFileAddress = catStrStr(rootAddress, "data/energyFile/");
const char *tempFileAddress   = catStrStr(energyFileAddress, "tempFile/");
const char *defaultFileAddress = catStrStr(energyFileAddress, "defaultFile/");
const char *QUACKoutFileAddress= catStrStr(energyFileAddress, "QUACKoutFile/");
const char *calFileAddress    = catStrStr(energyFileAddress, "calFile/");
const char *charmmFileAddress = catStrStr(rootAddress, "charmm/exec/gnu/");
const char *refine_1Address   = catStrStr(rootAddress, "refine/");
const char *TM_scoreAddress   = catStrStr(rootAddress, "TM_score/");
const char *scoreAddress      = catStrStr(rootAddress, "rosetta/main/source/bin/score.linuxgccrelease");
const char *databaseAddress   = catStrStr(rootAddress, "rosetta/main/database/");
const char *calRWplusAddress  = catStrStr(rootAddress, "calRWplus/");
const char *mybinAddress      = catStrStr(rootAddress, "mybin/");
const char *strideAddress     = catStrStr(rootAddress, "stride/");
const char *TM_alignAddress = catStrStr(rootAddress, "TMalign/");

//const char *answerAddress     = catStrStr(rootAddress, "data/answer/newAnswer/");
//const char *charmmFileAddress = catStrStr(energyFileAddress, "charmmFile/");
//const char inputAddress[]="/home/ws/zzZyj/MOPSO/data/input/test_thread";
//const char *inputAddress      = catStrStr(rootAddress, "data/input/originpdb/");
//const char *mybinAddress      = catStrStr(rootAddress, "mybin/");

// MOPSO Settings
const int tidSize =8;
const int objectiveNumber = 3; //  Multiple Obejectives settings, the size of objective function
const double TM_scoreThreshold = 0.13;
const int nRep = 50;                // Repository Size
const double Criterion = 0.000001;
const int lambdaLoopTimes = 1000;

// Problem Definition
const double angleMax = 180;
const double angleMin = -180;

const double phi1 = 2.05;
const double phi2 = 2.05;
const double phi = phi1 + phi2;
const double chi = 2 / ( phi - 2 + sqrt( phi * phi - 4 * phi ) );    // 0.73
const int bufferLen = 2048;

const double wMin = chi;                        // =chi  Inertia Weight
const double wMax = 1.2;
//const double c1 = chi*phi1;                 //  Personal Learning Coefficient
//const double c2 = chi*phi2;                 //  Global Learning Coefficient
const double c1max = 2.0;
const double c1min = 0.5;
const double c2max = 2.5;
const double c2min = 1.0;

const double PI = M_PI;
const double INF=100000000;
const double realMax = 100000000;
const double tiny = 0.0000000001;
