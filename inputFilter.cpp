//
// Created by advanced on 18-3-29.
//

#include "inputFilter.h"
#include "MOPSO.h"


int CountLines(char *filename) {
    ifstream ReadFile;
    int n=0;
    string tmp;
    ReadFile.open(filename,ios::in);//ios::in 表示以只读的方式读取文件
    if(ReadFile.fail())//文件打开失败:返回0
    {
        return 0;
    }
    else//文件存在
    {
        while(getline(ReadFile,tmp,'\n'))
        {
            n++;
        }
        ReadFile.close();
        return n;
    }
}

string ReadLine(char *filename,int line) {
    int lines,i=0;
    string temp;
    fstream file;
    file.open(filename,ios::in);
    lines=CountLines(filename);

    if(line<=0)
    {
        return "Error 1: 行数错误，不能为0或负数。";
    }
    if(file.fail())
    {
        return "Error 2: 文件不存在。";
    }
    if(line>lines)
    {
        return "Error 3: 行数超出文件长度。";
    }
    while(getline(file,temp)&&i<line-1)
    {
        i++;
    }
    file.close();
    return temp;
}

float mean(float a[],int b) {
    float c=0;
    for(int i=0;i<b;i++) {
        c = a[i] + c;
    }
    float mean=c/b;
    return mean;
}

double deviations(float x[],int y) {
    float c=0,k=0;
    for(int i=0;i<y;i++) {
        c = x[i] + c;
    }
    float mean=c/y;
    for(int i=0;i<y;i++){
        k=(x[i]-mean)*(x[i]-mean)+k;
    }
    double deviations=sqrt(k/y);
    return deviations;
}

void get_RSMD_score(const int &repID, const int &producedRepAmount) {
    //cout << "repID : " << repID << endl;

    char *RSMDoutput = catStrStr(answerAddress, "RSMDoriginFile.txt");                 //"/home/advanced/TR829_1.txt";
    for(int i = 0 ; i < producedRepAmount; i++) {
        const char *TMscoreAddress = catStrStr(rootAddress, "TMscore ");
        const char *multiPDB_address = catStrIntStr(answerAddress, "originCopy", i, ".pdb ");
        const char *originPDB_address = catStrIntStr(inputAddress, repID, ".pdb");

        const char *command = catStrStr(TMscoreAddress, multiPDB_address, originPDB_address, " >>", RSMDoutput);
       // cout << "command : " << command << endl << endl;
        //sprintf(TMscore,
          //      "  /home/advanced/实验环境/mopso/TMscore /home/advanced/TR%dAnswer/rep%d.pdb /home/advanced/TR%d/1.pdb >>/home/advanced/TR%d_1.txt",
            //    k,i,k,k);

        system(command);
    }

    int lines;
    string m;
    string str1,str2,str3,str4,str5,str6;
    lines=CountLines(RSMDoutput);
    for(int i=15;i<=lines;i=i+34)
    {
        m=ReadLine(RSMDoutput,i);
        istringstream is(m);//字符串分割
        is>>str1>>str2>>str3>>str4>>str5>>str6;
        //cout<<str3<<endl;
        ofstream outfile;
        char *finalRSMDfile = catStrStr(answerAddress, "RSMDfile.txt");

        outfile.open(finalRSMDfile,ios::app);
        outfile<<str6<<endl;
        outfile.close();

    }

    removeFile(RSMDoutput);
}

void get_RWPlus_score(const int &multiAmount) {
    char *RWplusOutput = catStrStr(answerAddress, "RWplusOriginFile.txt");                 //"/home/advanced/TR829_1.txt";
    //cout << "I'm here \n";

    pthread_t tids[tidSize];
    RW_Data *rw_data = new RW_Data[multiAmount];
    for (int i = 0; i < multiAmount; i++){
        rw_data[i].order = i;
    }

    int haveRun=0;
    while (haveRun < multiAmount){
        for (int i = 0; i < tidSize && haveRun + i < multiAmount; i++){
            pthread_create(&tids[i], NULL, getA_RwPlus_score, (void *) &rw_data[i]);
        }
        for (int i = 0; i < tidSize && haveRun + i < multiAmount; i++){
            pthread_join(tids[i], NULL);
        }
        haveRun += tidSize;
    }
    delete []rw_data;


    int lines;
    string m;
    string str1,str2,str3,str4,str5;
    lines=CountLines(RWplusOutput);
    for(int i=1;i<=lines;i=i+1)
    {
        m=ReadLine(RWplusOutput, i);
        istringstream is(m);//字符串分割
        is>>str1>>str2>>str3>>str4>>str5;
        //cout<<str3<<endl;
        ofstream outfile;
        char *RWplusFinal = catStrStr(answerAddress, "RWplusFile.txt");
        outfile.open(RWplusFinal, ios::app);
        outfile<<str4<<endl;
        outfile.close();

    }

    removeFile(RWplusOutput);
}

void *getA_RwPlus_score(void *p){
    RW_Data *rw_Data = (RW_Data *) p;
    char *RWplusOutput = catStrStr(answerAddress, "RWplusOriginFile.txt");                 //"/home/advanced/TR829_1.txt";
    const char *RWplusAddress = catStrStr(rootAddress, "calRWplus/");
    const char *multiPDB_address = catStrIntStr(answerAddress, "originCopy", rw_Data->order, ".pdb ");
    const char *command = catStrStr("cd ", RWplusAddress," && ./calRWplus ", multiPDB_address, " >>", RWplusOutput);
    //cout << "command : " << command << endl << endl;
    //sprintf(RWplus,
    //      "  /home/advanced/实验环境/mopso/calRWplus/calRWplus /home/advanced/TR%dAnswer/rep%d.pdb >>/home/advanced/TR%d.txt",
    //    k, i, k);

    system(command);
}

int *inputFilter(const int &repID, const int &multiAmount, const int &needAmount){
    int *sievedRes = new int[multiAmount];
    int *res = new int [needAmount];
    

    get_RSMD_score(repID, multiAmount);
    get_RWPlus_score(multiAmount);

    float *RWplus, *RSMD;
    RWplus = new float[multiAmount];
    RSMD = new float[multiAmount];

    ifstream infile;
    char *RWplusFinal = catStrStr(answerAddress, "RWplusFile.txt");
    infile.open(RWplusFinal);
    if(!infile) cout<<"error"<<endl;
    float t1;
    cout<<"存入数组RWplus"<<endl;
    float*p=&RWplus[0];
    while(infile>>t1)             //遇到空白符结束
    {
        *p=t1;
        p++;
    }
    infile.close();

    ifstream infile2;
    char *finalRSMDfile = catStrStr(answerAddress, "RSMDfile.txt");
    infile2.open(finalRSMDfile);
    if(!infile2) cout<<"error"<<endl;
    float t2;
    cout<<"存入数组RSMD1"<<endl;
    float*q=&RSMD[0];
    while(infile2>>t2)             //遇到空白符结束
    {
        *q=t2;
        q++;
    }
    infile2.close();

    double *rwplus,*rsmd;
    rwplus = new double[multiAmount];
    rsmd   = new double[multiAmount];
    for (int i = 0; i < multiAmount; i++)
    {
        rwplus[i] = RWplus[i] - mean(RWplus,multiAmount);
        rwplus[i] = rwplus[i] / deviations(RWplus,multiAmount);
        rsmd[i] = RSMD[i] - mean(RSMD,multiAmount);
        rsmd[i] = rsmd[i] / deviations(RSMD,multiAmount);
         cout<<i<<"  "<<RSMD[i]<<endl;
         //cout<<i<<"  "<<rwplus[i]<<endl;
    }
    double row = 1, gama = 60;
    int inc = 0;
    for (int i = 0 ; i < multiAmount; i++)
    {
        if(rsmd[i]*rsmd[i]+rwplus[i]*rwplus[i]>=row&&(acos((-0.5*rwplus[i]-sqrt(0.75)*rsmd[i])/sqrt(rsmd[i]*rsmd[i]+rwplus[i]*rwplus[i])))*180/3.1415926<=gama)
        {
            sievedRes[inc++] = i;
        }
    }

    if (inc < needAmount){
        cout << "Not satisfy need amount. rep : " << repID << endl;
        exit(-1);
    }

    for (int i = 0; i < needAmount; i++)
        res[i] = sievedRes[i];

    cout << "The result is : \n";
    for (int i = 0; i < needAmount; i++)
        cout << res[i] << " ";
    cout << endl;

    for (int j = 0; j < multiAmount; j++){
        //const char *multiPDB_address = catStrIntStr(answerAddress, "originCopy", 1, ".pdb ");
        //removeFile(multiPDB_address);
    }
    removeFile(RWplusFinal);
    removeFile(finalRSMDfile);
    return sievedRes;
}


