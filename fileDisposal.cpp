//
// Created by ws on 16-7-15.
//

#include "MOPSO.h"



int getLines(char *address){
    // get the lines size of a file
    ifstream inFile(address);
    if (!inFile.is_open()){
        cout << "can't open " << address << endl;
        exit(0);
    }
    char buffer[bufferLen];
    int cntLines = 0;
    while (! inFile.eof()) { inFile.getline(buffer, bufferLen); cntLines++; }
    cntLines--;
    inFile.close();
    return cntLines;
}

void getAtom(Atom *atom, int cntLines, int &atomNum, char *address) {
    freopen(address, "r", stdin);
    char buffer[bufferLen];
    char judge[bufferLen];
    atomNum = 0;
    for (int i=0; i<cntLines; i++){
        strcpy(judge, "");              //  bug : if the last line is "", then judge is not read a new "", but remain the old value
        scanf("%s", judge);

        if (strcmp(judge, "ATOM") == 0){
            scanf("%d%s%s%d%lf%lf%lf",
                  &atom[atomNum].firNumber,atom[atomNum].name, atom[atomNum].group,
                  &atom[atomNum].secNumber, &atom[atomNum].x, &atom[atomNum].y, &atom[atomNum].z);
            //for(char c=getchar();c!='\n';c=getchar()) str[++n]=c;
            scanf("%[^\n]", buffer);
            atomNum++;
        }
        else
            scanf("%[^\n]", buffer);
    }
}

char * inputSeq(){
    // read in the sequence
    char *seq;
    char buffer[bufferLen];
    char *filename;
    ifstream inFile;
    filename = catStrStr(inputAddress,"seq.txt");
    openFile(filename, inFile);
    int cntLines = 0;
    while (! inFile.eof()) { inFile.getline(buffer, bufferLen); cntLines++; }
    cntLines--;
    inFile.close();

    openFile(filename, inFile);
    char tmp[cntLines][bufferLen];
    int len=0;
    inFile.getline(buffer, bufferLen);                    //  the first line is not needed
    for (int i=0; i<cntLines-1; i++) {
        inFile.getline(tmp[i], bufferLen);
        len += strlen(tmp[i]);
    }

    seq = new char [len+1];
    strcpy(seq, "");
    for (int i=0; i<cntLines-1; i++)  strcat(seq, tmp[i]);
    inFile.close();
    return seq;
}

// input the coordinaries and angles

void openFile(const char fileDir[10000], ifstream &infile){
    infile.open(fileDir);
    if (!infile) {
        cout << "can't open:::" << fileDir << endl;
        exit(-1);
    }
}

void removeFile(const char * str){
    int k = remove(str);
    if (k==-1) {
        cout << "Fail to delete:::  " << str << endl;
        exit(-1);
    }
}

void removeMultiplyFile_afterInput(){
    char *seqAddress = catStrStr(inputAddress, "seq.txt");
    removeFile(seqAddress);
    for (int i = 0; i < inputSize; i++){
        char *prePhiAddress = catStrIntStr(inputAddress, "prePhi", i+1, ".txt");
        removeFile(prePhiAddress);
    }
    for (int i = 0; i < Population; i++){
        char *particleAddress = catStrIntStr(inputAddress, "particleO_", i+1, ".txt");
        char *phiAddress = catStrIntStr(inputAddress, "phi", i+1, ".txt");
        removeFile(particleAddress);
        removeFile(phiAddress);
    }

}

void printParticleCost(Particle * particle, int iterator){
    char *costFile = catStrStr(answerAddress, "cost.txt");
    FILE *outFile;
    outFile = fopen(costFile, "a");
    if (outFile == NULL) {
        cout << "Can't print rep answer" << endl;
        exit(-1);
    }
    fprintf(outFile, "iterator=%d\n", iterator);
    for (int i=0; i<Population; i++){
        fprintf(outFile, "i=%d ", i);
        for (int j=0; j<objectiveNumber; j++){
            fprintf(outFile, "cost%d =%10.3f    ", j, particle[i].Cost[j]);
        }
        fprintf(outFile, "\n");
    }
    fclose(outFile);
}

//  fprintf
void outputAnswer(myRep &rep) {
    //sortRepByAngle(rep, angle);
    sortRepByLambda(rep);
    myRep ::iterator it1;
    int repSeqNumber=0;
    int size = static_cast<int>(rep.size());
    int loop = size;
    for (int q=0; q < loop; q++){
        it1 = rep.begin();
        int cnt = sortAns[q];
        for (int p = 0; p < cnt; p++) {
            it1++;
        }
        printPdb(it1, repSeqNumber);
        char *filename;
        filename = catStrIntStr(answerAddress, "rep", repSeqNumber, ".txt");
        FILE *outFile;
        outFile = fopen(filename, "w");
        if (outFile == NULL) {
            cout << "Can't print rep answer" << endl;
            exit(-1);
        }
        fprintf(outFile, "iterator=%d\n",it1->iterator);
        for (int i = 0; i < objectiveNumber; i++){
            fprintf(outFile, "cost%d =%10.3f ", i, it1->Cost[i]);
        }
        fprintf(outFile, "\n");
        int len = it1->sizeOfAddO_origin;
        int t=0;
        for (int j=0; j<len; j+=3){
            fprintf(outFile, "%4d ", t);
            switch (t % 4){
                case 0: fprintf(outFile, " N  "); break;
                case 1: fprintf(outFile, " CA "); break;
                case 2: fprintf(outFile, " C  "); break;
                case 3: fprintf(outFile, " O  "); break;
                default:
                    break;
            }
            fprintf(outFile, "%8.3f %8.3f %8.3f\n",
                    it1->addO_origin[j], it1->addO_origin[j+1], it1->addO_origin[j+2]);
            t++;
        }
        fclose(outFile);

        // output all rep's cost to a single file
        filename = catStrStr(answerAddress, "repCost.txt");
        outFile = fopen(filename, "a");
        fprintf(outFile, "%d ", repSeqNumber);
        for (int i = 0; i < objectiveNumber; i++){
            fprintf(outFile, "%10.3f ", it1->Cost[i]);
        }
        fprintf(outFile, "\n");
        fclose(outFile);

        repSeqNumber++;
    }
    //add the angle_select 
    
}



void printTime(const int choice, const int loopTimes) {
    getEndTime();
    const char *timeFile = catStrStr(answerAddress, "time.txt");
    ofstream outfile(timeFile, ios::app);
    long duration;
    duration = endTime-starTime;
    switch (choice){
        case 0: outfile << "TOTAL Time: \n"; break;
        case 1: outfile << "Have Run " << loopTimes << " times : "; break;
        case 2: outfile << "Cost 2 : "; break;
        case 3: outfile << "Transfer Time: "; break;
        case 4: outfile << "Input Time: "; break;
        case 5: outfile << "Adaption Time: "; break;
        case 6: outfile << "IntoRep Time: "; break;
        default: break;
    }
    long seconds = duration % 60;
    duration /= 60;
    long minutes = duration % 60;
    duration /= 60;
    long hours = duration;

    outfile << hours << ":" << minutes << ":" << seconds << endl;
    outfile.close();
}

void printCostTime(long rosettaTime, long calTime , long charmmTime, long mybinTime){
    const char *timeFile = catStrStr(answerAddress, "EnergyTime.txt");
    ofstream outfile(timeFile);

    outfile << "rosetta : " << rosettaTime << endl;
    outfile << "cal     : " << calTime << endl;
    outfile << "charmm  : " << charmmTime << endl;
    outfile << "mybin   : " << mybinTime << endl;

    outfile.close();

}

void createNewFold(){
    const char *checkExist = catStrStr("if [ -d ", answerAddress, " ]; then\nrm -r ",
                                       answerAddress, "\nfi");
    system(checkExist);
    const char *mkdirCommand = catStrStr("mkdir ",answerAddress);
    int t= system(mkdirCommand);
    if (t != 0){
        cout << "can't create a new fold for answers" << endl;
        exit(-1);
    }

}

void printPdb(Particle * particle){
    FILE * output;
    const char *fileName = catStrIntStr(tempFileAddress, "temp", particle->index, ".pdb");
    output = fopen(fileName, "w");
    if (output == NULL) {
        cout << "Can't write temp.pdb" << endl;
        exit(-1);
    }
    int numID = 0;
    int lenOfAA = particle->numAA;
    for (int i = 0; i < lenOfAA; i++){
        char *AA = getRelativeAA(particle->seq[i]);
        int k1 = i*4*3, k2 = i*4*3+3, k3 = i*4*3+6, k4 = i*4*3+9;
        fprintf(output, "ATOM   %5d N   %3s  %4d    %8.3f%8.3f%8.3f  1.00  0.00\n",
        //        par_firNumber[numID],   AA, par_secNumber[numID  ], particle->addO_origin[k1], particle->addO_origin[k1+1], particle->addO_origin[k1+2]);
                4*i+1, AA, i+1, particle->addO_origin[k1], particle->addO_origin[k1+1], particle->addO_origin[k1+2]);
        fprintf(output, "ATOM   %5d CA  %3s  %4d    %8.3f%8.3f%8.3f  1.00  0.00\n",
                //par_firNumber[numID+1], AA, par_secNumber[numID+1], particle->addO_origin[k2], particle->addO_origin[k2+1], particle->addO_origin[k2+2]);
                4*i+2, AA, i+1, particle->addO_origin[k2], particle->addO_origin[k2+1], particle->addO_origin[k2+2]);
        fprintf(output, "ATOM   %5d C   %3s  %4d    %8.3f%8.3f%8.3f  1.00  0.00\n",
                //par_firNumber[numID+2], AA, par_secNumber[numID+2], particle->addO_origin[k3], particle->addO_origin[k3+1], particle->addO_origin[k3+2]);
                4*i+3, AA, i+1, particle->addO_origin[k3], particle->addO_origin[k3+1], particle->addO_origin[k3+2]);
        fprintf(output, "ATOM   %5d O   %3s  %4d    %8.3f%8.3f%8.3f  1.00  0.00\n",
                //par_firNumber[numID+3], AA, par_secNumber[numID+3], particle->addO_origin[k4], particle->addO_origin[k4+1], particle->addO_origin[k4+2]);
                4*i+4, AA, i+1, particle->addO_origin[k4], particle->addO_origin[k4+1], particle->addO_origin[k4+2]);
        delete AA;
        numID+=4;
    }
    fclose(output);
}

void printPdb(Particle * particle, const char *fileName){
    FILE * output;
    output = fopen(fileName, "w");
    if (output == NULL) {
        cout << "Can't write pdb : "<< fileName << endl;
        exit(-1);
    }
    int lenOfAA = particle->numAA;
    int numID = 0;
    for (int i = 0; i < lenOfAA; i++){
        char *AA = getRelativeAA(particle->seq[i]);
        int k1 = i*4*3, k2 = i*4*3+3, k3 = i*4*3+6, k4 = i*4*3+9;
        fprintf(output, "ATOM   %5d N   %3s  %4d    %8.3f%8.3f%8.3f  1.00  0.00\n",
                par_firNumber[numID],   AA, par_secNumber[numID  ], particle->addO_origin[k1], particle->addO_origin[k1+1], particle->addO_origin[k1+2]);
        fprintf(output, "ATOM   %5d CA  %3s  %4d    %8.3f%8.3f%8.3f  1.00  0.00\n",
                par_firNumber[numID+1], AA, par_secNumber[numID+1], particle->addO_origin[k2], particle->addO_origin[k2+1], particle->addO_origin[k2+2]);
        fprintf(output, "ATOM   %5d C   %3s  %4d    %8.3f%8.3f%8.3f  1.00  0.00\n",
                par_firNumber[numID+2], AA, par_secNumber[numID+2], particle->addO_origin[k3], particle->addO_origin[k3+1], particle->addO_origin[k3+2]);
        fprintf(output, "ATOM   %5d O   %3s  %4d    %8.3f%8.3f%8.3f  1.00  0.00\n",
                par_firNumber[numID+3], AA, par_secNumber[numID+3], particle->addO_origin[k4], particle->addO_origin[k4+1], particle->addO_origin[k4+2]);
        delete AA;
        numID+=4;
    }
    fclose(output);
}

void printCopiedOrigin(const Particle &particle, const int repNum) {
    // generate particleO_21.txt ...
    char *filename;
    filename = catStrIntStr(inputAddress, "particleO_", repNum+1, ".txt");

    ofstream outFile(filename);
    int len = particle.sizeOfAddO_origin;
    for (int i=0; i<len; i++){
        outFile << particle.addO_origin[i] << " ";
        if (i % 3 == 2) outFile << endl;
    }

    // generate phi21.txt ...
    filename = catStrIntStr(inputAddress, "phi", repNum+1, ".txt");
    //outFile.open(filename);
    ofstream outFileForPhi(filename);
    len = 3*particle.numAA-5;
    for (int j=0; j < len; j+=3){          // last : j=3*numAA-6; j+2=3*numAA-4
        outFileForPhi << particle.Position[j] << " ";
        outFileForPhi << particle.Position[j+2] << endl;
    }
}

void printPdb(list<Rep>::iterator it, const int &repNum) {
    //cout << "seq : " << it->seq << endl;

    FILE * output;
    const char *filename = catStrIntStr(answerAddress, "rep", repNum, ".pdb");
    output = fopen(filename, "w");
    if (output == NULL) {
        cout << "Can't write temp.pdb" << endl;
        exit(-1);
    }
    int lenOfAA = it->numAA;
    int numID = 0;
    int AA_ID = 0;
    for (int i = 0; i < lenOfAA; i++){
        char *AA = getRelativeAA(it->seq[i]);
        int k1 = i*4*3, k2 = i*4*3+3, k3 = i*4*3+6, k4 = i*4*3+9;
        fprintf(output, "ATOM   %5d N   %3s  %4d    %8.3f%8.3f%8.3f  1.00  0.00\n",
                par_firNumber[numID]  , AA, par_secNumber[numID], it->addO_origin[k1], it->addO_origin[k1+1], it->addO_origin[k1+2]);
                //4*i+1, AA, AA_ID, it->addO_origin[k1], it->addO_origin[k1+1], it->addO_origin[k1+2]);
        fprintf(output, "ATOM   %5d CA  %3s  %4d    %8.3f%8.3f%8.3f  1.00  0.00\n",
                par_firNumber[numID+1], AA, par_secNumber[numID+1], it->addO_origin[k2], it->addO_origin[k2+1], it->addO_origin[k2+2]);
                //4*i+2, AA, AA_ID, it->addO_origin[k2], it->addO_origin[k2+1], it->addO_origin[k2+2]);
        fprintf(output, "ATOM   %5d C   %3s  %4d    %8.3f%8.3f%8.3f  1.00  0.00\n",
                par_firNumber[numID+2], AA, par_secNumber[numID+2], it->addO_origin[k3], it->addO_origin[k3+1], it->addO_origin[k3+2]);
                //4*i+3, AA, AA_ID, it->addO_origin[k3], it->addO_origin[k3+1], it->addO_origin[k3+2]);
        fprintf(output, "ATOM   %5d O   %3s  %4d    %8.3f%8.3f%8.3f  1.00  0.00\n",
                par_firNumber[numID+3], AA, par_secNumber[numID+3], it->addO_origin[k4], it->addO_origin[k4+1], it->addO_origin[k4+2]);
                //4*i+4, AA, AA_ID, it->addO_origin[k4], it->addO_origin[k4+1], it->addO_origin[k4+2]);
        delete AA;
        AA_ID++;
        numID += 4;
    }
    fclose(output);
}

void runScwrl(int index){
    const char *fileName = catStrIntStr(tempFileAddress, "temp", index, ".pdb");
    const char *newFile = catStrIntStr(tempFileAddress,"temp_", index, ".pdb");
    const char *command = catStrStr(refine_1Address, "scwrl4/Scwrl4 -i ", fileName,
                                    " -o ", newFile);
    system(command);
}

void removeTmp(int index){
    const char *fileName = catStrIntStr(tempFileAddress, "temp", index, ".pdb");
    const char *newFile = catStrIntStr(tempFileAddress,"temp_", index, ".pdb");
    removeFile(fileName);
    removeFile(newFile);
}

double getCost_rosetta(int i) {
    //remove("/home/ws/zzZyj/data/temp1.pdb"); remove("default.sc");
    //printf("Particle ID = %d, Enter getCost_rosetta\n", i);
    const char *temp_File = catStrIntStr(tempFileAddress, "temp_", i, ".pdb");
    const char *defaultFile = catStrIntStr(defaultFileAddress, "default", i, ".sc");
    const char *command = catStrStr(scoreAddress, " -database ", databaseAddress," -s ",
                                    temp_File, " -out:file:scorefile ", defaultFile);
    system(command);
    ifstream infile;
    openFile(defaultFile, infile);
    char buffer[10000], strAns[10000];
    double ans=0;
    infile.getline(buffer, 10000);
    infile >> buffer >> strAns;
    int len = static_cast<int> (strlen(strAns));
    if (strAns[0]=='n' || (len>=2 && strAns[1]=='n'))
        ans = INF;
    else
        ans = atof(strAns);
    infile.close();
    //cout << "defaultFile = " << defaultFile << endl;
    removeFile(defaultFile);
    //printf("Particle ID = %d, Out getCost_rosetta\n", i);
    return ans;
}

double getCost_calRWplus(int i){
    //printf("Particle ID = %d, Enter getCost_calRWplus\n", i);
    int r_index,j,k;
 //char *tempFileAddress   = "/home/gengling/MOPSO/data/energyFile/tempFile/";
 //char *QUACKoutFileAddress=  "/home/gengling/MOPSO/data/energyFile/QUACKoutFile/";
    double ans=0;
    while (ans==0){
        char buffer[bufferLen],cost[20];
        const char *str = catStrStr("cd ", calRWplusAddress, " && ./calRWplus ",tempFileAddress,"temp_");
        const char *command1 = catStrIntStr(str, i, ".pdb > ");
        const char *command2 = catStrIntStr(calFileAddress, "cal", i, ".txt");
        const char *command = catStrStr(command1, command2);
        system(command);
        ifstream infile;
        openFile(command2, infile);
        infile.getline(buffer, bufferLen);
        for(j=0;j<strlen(buffer);j++)
        {
            if(buffer[j]=='k') r_index=j;
        }
         k=0;
         for(j=14;j<r_index;j++)
        {
            cost[k]=buffer[j];k++;
        }
        cost[k]='\0';
        sscanf(cost,"%lf",&ans);
        infile.close();
        removeFile(command2);
        delete command;  delete command1; delete command2;
    }
    //printf("Particle ID = %d, Out getCost_calRWplus\n", i);
    return ans;
}


double getCost_charmm(int i){
    //printf("Particle ID = %d, Enter getCost_charmm\n", i);
    const char *tmpI = catStrIntStr("tmp", i, "");
    const char *temp_File = catStrIntStr(tempFileAddress, "temp_", i, ".pdb");
    const char *charmmFile = charmmFileAddress;
    const char *awkFile = catStrStr(charmmFile, "fixpdb.awk");
    const char *charmmTempFile = catStrStr(charmmFile, tmpI, ".pdb");
    const char *awkCommand = catStrStr("awk -f ", awkFile, " segid=PROT <",
                                       temp_File, "> ", charmmTempFile);
    system(awkCommand);
   // cout << "run awk"<< endl;
    // rewrite template.inp as tmp1.inp
    const char *tmpInp = catStrStr(charmmFile, tmpI, ".inp");
    const char *templateFile = catStrStr(charmmFile, "template.inp");
    rewriteTemplate(i, templateFile, tmpInp, tmpI);


    const char *tmpOut = catStrStr(charmmFile, tmpI, ".out");
    const char *getEnergyFile = catStrStr("cd ", charmmFileAddress," && ./charmm <",
                                        tmpInp, "> ", tmpOut);
    system(getEnergyFile);
    // read the value of energy in the energyFile
    double ans = getEnergyValue(tmpOut);

    // remove useless file
    removeFile(tmpInp);
    removeFile(tmpOut);
    //printf("Particle ID = %d, Out getCost_charmm\n", i);
    return ans;
}
/*
double getCost_mybin(int i){
    const char *command01 = catStrIntStr(tempFileAddress, "temp_", i, ".pdb");
    const char *command02 = catStrIntStr(mybinAddress, "temp_", i, ".pdb");
    const char *command0 = catStrStr("cp ", command01, " ", command02);
    system(command0);

    double ans=0;
    while (ans==0){
        int lineNum=0;
        char buffer[10000];
        //system("/home/ws/GL/mybin/calcquarkenergy /home/ws/GL/mybin/ /home/ws/GL/mybin/temp1.pdb>>QUACKout.txt");
        const char *str = catStrStr(mybinAddress, "calcquarkenergy ", mybinAddress, " ",
                                    mybinAddress, "temp_");
        const char *command1 = catStrIntStr(str, i, ".pdb>>");
        const char *command2 = catStrIntStr(QUACKoutFileAddress, "QUACKout", i, ".txt");
        const char *command = catStrStr(command1, command2);
        system(command);
        ifstream infile;
        openFile(command2, infile);
        while (infile.getline(buffer, 10000)) lineNum++;
        infile.close();
        openFile(command2, infile);
        for (int i=0; i<lineNum-1; i++) infile.getline(buffer, 10000);
        infile >> ans;
        infile.close();
        removeFile(command2);
        delete command;  delete command1; delete command2;
    }
    removeFile(command02);

    delete command0;  delete command01; delete command02;

    return ans;
}

*/
void rewriteTemplate(int i, const char *templateFile, const char *tmpIFile, const char *tmpI) {
    ifstream infile(templateFile);
    ofstream outfile(tmpIFile, ios::out);
    char buffer[bufferLen];
    int cntLines=0;
    while (!infile.eof()){
        infile.getline(buffer, bufferLen);
        cntLines++;
        if (cntLines == 21 || cntLines==33){
            strReplace(buffer, tmpI, 27, 8);
        }
        outfile << buffer << endl;

    }
    infile.close();
    outfile.close();
}

double getEnergyValue(const char *fileName){
    ifstream energyFile(fileName);
    char buffer[bufferLen];
    double ans=0;
    while (1){
        energyFile >> buffer;
        if (strcmp(buffer, "ENER>") == 0){
            energyFile >> buffer;
            energyFile >> buffer;
            ans = atof(buffer);
            break;
        }
    }
    return ans;
}

void printCurrentId(int id){
    cout << "current ID : " << id << endl;
}

void printParticle(Particle * particle){
    FILE * output;
    const char *fileName = catStrIntStr("/home/advanced/MOPSO/data/input/", "MOPSO_Particle", particle->index, ".pdb");
    output = fopen(fileName, "w");
    if (output == NULL) {
        cout << "Can't write temp.pdb" << endl;
        exit(-1);
    }
    int lenOfAA = particle->numAA;
    for (int i = 0; i < lenOfAA; i++){
        char *AA = getRelativeAA(particle->seq[i]);
        int k1 = i*4*3, k2 = i*4*3+3, k3 = i*4*3+6, k4 = i*4*3+9;
        fprintf(output, "ATOM   %4d  N   %3s   %3d    %8.3f%8.3f%8.3f  1.00  0.00\n",
                4*i+1, AA, i+1, particle->addO_origin[k1], particle->addO_origin[k1+1], particle->addO_origin[k1+2]);
        fprintf(output, "ATOM   %4d  CA  %3s   %3d    %8.3f%8.3f%8.3f  1.00  0.00\n",
                4*i+2, AA, i+1, particle->addO_origin[k2], particle->addO_origin[k2+1], particle->addO_origin[k2+2]);
        fprintf(output, "ATOM   %4d  C   %3s   %3d    %8.3f%8.3f%8.3f  1.00  0.00\n",
                4*i+3, AA, i+1, particle->addO_origin[k3], particle->addO_origin[k3+1], particle->addO_origin[k3+2]);
        fprintf(output, "ATOM   %4d  O   %3s   %3d    %8.3f%8.3f%8.3f  1.00  0.00\n",
                4*i+4, AA, i+1, particle->addO_origin[k4], particle->addO_origin[k4+1], particle->addO_origin[k4+2]);
        delete AA;
    }
    fclose(output);
}
