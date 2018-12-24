//
// Created by advanced on 17-1-16.



/*
void printRep(myRep & rep, int len){
    int repNum=0;
    repNum++;
    char *filename;
    filename = catStrIntStr("/home/gengling/MOPSO/data/rep/rep", repNum, ".txt");
    ofstream outFile(filename);
    myRep ::iterator it2;
    int rr=0;
    for (it2=rep.begin(); it2!= rep.end(); it2++){
        rr++;
        outFile <<  rr << " cost[0]=" << it2->Cost[0] << " cost[1]=" << it2->Cost[1] << endl;
        outFile << "Following is the x,y,z coordinary: "<< endl;
        for (int i=0; i<len; i+=3){
            outFile<< i/3+1 << " " << it2->addO_origin[i] << " " << it2->addO_origin[i+1] << " " << it2->addO_origin[i+2]<< endl;
        }
        outFile << "------------------------------------------------------------------"<< endl;
    }
    delete filename;
    outFile.close();
}
 */