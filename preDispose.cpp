//
// Created by advanced on 17-1-16.
//

#include "MOPSO.h"

int inputSize;               // Population = inputFilesNumber + multiplyFilesNumber
char *inputAddress;
char *answerAddress;
char *logAddress;
char *energyOutputAddress;
int  MaxIt;
int multiplyNumber;
char **argv;
int Population;
int startNum;
int numAA;
int currentIteration_Times;
int times_for_each_play;
double VEL_SMALL_RANGE;
double VEL_BIG_RANGE;

void preDisposeInputParametersAndFiles(char **argv) {
    getInputParameter(argv);
    inputAddress = catStrStr(argv[1], "/");
    answerAddress = catStrStr(argv[1], "Answer/");
    logAddress = catStrStr(answerAddress, "log.txt");
    energyOutputAddress = catStrStr(answerAddress, "energyOutput.txt");

    createNewFold();
    disposePDB();
    multiplyNumber = Population - inputSize;
    getParameters_for_TMalign();
}

void disposePDB(){
    createSeqTxt();
    for (int i=0;  i < inputSize; i++){
        createParticleTxt(i+1);
        createPhi(i+1);
    }


}

void createSeqTxt(){
    //char address[100] = "/home/ws/zzZyj/MOPSO/data/temp.pdb";
    char *address = catStrIntStr(inputAddress, 1, ".pdb");
    int cntLines = getLines(address);

    // input the data in pdb
    Atom *atom;
    int atomNum;
    atom = new Atom[cntLines+100];
    getAtom(atom, cntLines, atomNum, address);

    // create seq.txt
    char *seqAddress = catStrStr(inputAddress, "seq.txt");
    freopen(seqAddress, "w", stdout);
    printf(">TR233\n");
    int k = 0;
    for (int i=0; i<atomNum; i++){
        if (atom[i].secNumber != k){
            printf("%c", getAbbreviation(atom[i].group));
            k = atom[i].secNumber;
        }
    }
    printf("\n");
    freopen(logAddress, "a", stdout);

}

void createParticleTxt(int k){
    char *address = catStrIntStr(inputAddress, k, ".pdb");
    int cntLines = getLines(address);

    // input the data in pdb
    Atom *atom;
    atom = new Atom[cntLines+100];
    int atomNum;
    getAtom(atom, cntLines, atomNum, address);

    //  creat particle.txt
    char *particleAddress = catStrIntStr(inputAddress, "particleO_", k, ".txt");
    freopen(particleAddress, "w", stdout);
    for (int i = 0; i < atomNum; i++){
        if (isImportantAtom(atom[i]))
            printf("%lf %lf %lf\n", atom[i].x-atom[0].x, atom[i].y-atom[0].y, atom[i].z-atom[0].z);
    }

    // print numbers in order to be the same with input pdb
    if (k == 1){
        int number_Length = 0;
        for (int i = 0; i < atomNum; i++){
            if (isImportantAtom(atom[i]))
                number_Length++;
        }
        par_firNumber = new int [number_Length];
        par_secNumber = new int [number_Length];
        int numID = 0;
        for (int i = 0; i < atomNum; i++)
            if (isImportantAtom(atom[i])){
                par_firNumber[numID] = atom[i].firNumber;
                par_secNumber[numID] = atom[i].secNumber;
                numID++;
            }

    }

    freopen(logAddress, "a", stdout);

}

void createPhi(int k){
    //system(" cd /home/ws/GL/stride && ./stride -mFile /home/ws/zzZyj/MOPSO/data/temp.pdb > /home/ws/zzZyj/MOPSO/data/tempPhi.txt");
    char *pdbAddress = catStrIntStr(inputAddress, k, ".pdb");
    char *prePhiAddress = catStrIntStr(inputAddress, "prePhi", k, ".txt");
    char *runStride = catStrStr("cd ", strideAddress, " && ./stride -mFile ",
                                pdbAddress, " > ", prePhiAddress);
    //  stride.exe
    system(runStride);
    // getLines
    int cntLines = getLines(prePhiAddress);

    //  create Phi.txt
    freopen(prePhiAddress, "r", stdin);
    char *phiAddress = catStrIntStr(inputAddress, "phi", k, ".txt");
    freopen(phiAddress, "w", stdout);
    char firstWord[100];
    int lineIt = 1;
    scanf("%s", firstWord);
    while (strcmp(firstWord, "ASG") != 0){
        char *buffer;
        buffer = new char[1000];
        scanf("%[^\n]", buffer);
        lineIt++;
        scanf("%s", firstWord);
    }
    // cout phi (ignore the first 360 and the last 360)
    double data1, data2;
    // first one
    scanf("%*s%*s%*s%*s%*s%*s%lf%lf%*f%*s", &data1, &data2);
    printf("%lf\n", data2);
    lineIt++;
    // ordinary ones
    for (int i = lineIt; i<cntLines; i++){
        scanf("%*s%*s%*s%*s%*s%*s%*s%lf%lf%*f%*s", &data1, &data2);
        printf("%lf %lf\n", data1, data2);
    }
    // last one
    scanf("%*s%*s%*s%*s%*s%*s%*s%lf%lf%*f%*s", &data1, &data2);
    printf("%lf\n", data1);
    freopen(logAddress, "a", stdout);
}

bool isImportantAtom(const Atom atom){
    const char *name = atom.name;
    if (strcmp(name, "N")==0 || strcmp(name, "CA")==0 || strcmp(name, "C")==0 || strcmp(name, "O")==0)
        return 1;
    else
        return 0;
}

char getAbbreviation(char *str){
    if (strcmp(str, "ALA") == 0) return 'A';
    if (strcmp(str, "ARG") == 0) return 'R';
    if (strcmp(str, "ASN") == 0) return 'N';
    if (strcmp(str, "ASP") == 0) return 'D';
    if (strcmp(str, "CYS") == 0) return 'C';
    if (strcmp(str, "GLN") == 0) return 'Q';
    if (strcmp(str, "GLU") == 0) return 'E';
    if (strcmp(str, "GLY") == 0) return 'G';
    if (strcmp(str, "HIS") == 0) return 'H';
    if (strcmp(str, "ILE") == 0) return 'I';
    if (strcmp(str, "LEU") == 0) return 'L';
    if (strcmp(str, "LYS") == 0) return 'K';
    if (strcmp(str, "MET") == 0) return 'M';
    if (strcmp(str, "PHE") == 0) return 'F';
    if (strcmp(str, "PRO") == 0) return 'P';
    if (strcmp(str, "SER") == 0) return 'S';
    if (strcmp(str, "THR") == 0) return 'T';
    if (strcmp(str, "TRP") == 0) return 'W';
    if (strcmp(str, "TYR") == 0) return 'Y';
    if (strcmp(str, "VAL") == 0) return 'V';
}

void getArgv(){
    argv = new char *[100];
    for (int i=0; i<100; i++)
        argv[i] = new char [100];
    for (int i=1; i<4; i++) {
        cin >> argv[i];
    }
}

void getParameters_for_TMalign(){
    int tmp = getLines(catStrStr(inputAddress, "particleO_1.txt"));
    numAA = tmp/4;
    const int divideNumber[7] = {0, 0, 1, 3, 6, 10, 15};                // 6C2 =15
    int fileRange_for_TM_align = min(inputSize, 6);         // max for TM_align is 1.pdb, 2.pdb, ... , 6.pdb
    times_for_each_play = MaxIt / divideNumber[fileRange_for_TM_align];           // MaxIt = 3000, inputSize=7, fileRange = 6
    if (times_for_each_play == 0)        //  MaxIt < divideNumber
        times_for_each_play = 1;

    char *seq;
    seq = inputSeq();//input seq
    int len = static_cast<int>(strlen(seq));
    if (len > 250){
        VEL_BIG_RANGE = 0.8;
        VEL_SMALL_RANGE = 0.4;
    }
    else if (len > 120){
        VEL_BIG_RANGE = 1.0;
        VEL_SMALL_RANGE = 0.5;
    }
    else{
        VEL_BIG_RANGE = 1.2;
        VEL_SMALL_RANGE = 0.6;
    }
}
