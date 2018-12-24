#include "MOPSO.h"

char *catStrIntStr(const char *c1, int ss,const char *c2){
    int len1 = strlen(c1), len2 = strlen(c2);
    char *str3;
    str3 = new char [len1+len2+10];
    strcpy(str3, c1);
    char tmp[100];
    sprintf(tmp, "%d", ss);
    strcat(str3, tmp);
    strcat(str3, c2);

    return str3;
}

char *catStrIntStr(const char *c1, const char *c2, int ss, const char *c3){
    const char *str1 = catStrStr(c1, c2);
    return catStrIntStr(str1, ss, c3);
}

char *catStrStr(const char *c1, const char *c2){
    int len1 = strlen(c1), len2 = strlen(c2);
    char *str3;
    str3 = new char [len1+len2+10];
    strcpy(str3, c1);
    strcat(str3,c2);
    return str3;
}

char *catStrStr(const char *c1, const char *c2, const char *c3){
    char *str = catStrStr(c1, c2);
    return catStrStr(str, c3);
}

char *catStrStr(const char *c1, const char *c2, const char *c3, const char *c4){
    const char *str = catStrStr(c1, c2);
    return catStrStr(str, c3, c4);
}

char *catStrStr(const char *c1, const char *c2, const char *c3, const char *c4, const char *c5){
    const char *str = catStrStr(c1, c2);
    return catStrStr(str, c3, c4, c5);
}

char *catStrStr(CONCH *c1, CONCH *c2, CONCH *c3, CONCH *c4, CONCH *c5, CONCH *c6){
    const char *str = catStrStr(c1, c2);
    return catStrStr(str, c3, c4, c5, c6);
}

char *catStrStr(CONCH *c1, CONCH *c2, CONCH *c3, CONCH *c4, CONCH *c5, CONCH *c6, CONCH *c7){
    const char *str = catStrStr(c1, c2);
    return catStrStr(str, c3, c4, c5, c6, c7);
}

char *catStrStr(CONCH *c1, CONCH *c2, CONCH *c3, CONCH *c4, CONCH *c5, CONCH *c6, CONCH *c7, CONCH *c8){
    const char *str = catStrStr(c1, c2);
    return catStrStr(str, c3, c4, c5, c6, c7, c8);
}

void cpyDoubleArray(double *& p1, double *& p2, int n){                //  p1,p2 is Position,  copy p2 to p1
    for (int i=0; i<n; i++){
        p1[i]  = p2[i];

    }
}

void cpyIntArray(int *& p1, int *& p2, int n){
    for (int i=0; i<n; i++){
        p1[i]  = p2[i];

    }
}

void checkV(double *tmp, int n){
    for (int i=0; i<n; i++) {
        if (tmp[i] > VelMax[i]) tmp[i] = VelMax[i];
        else if (tmp[i] < -VelMax[i]) tmp[i] = -VelMax[i];
    }
}

void checkP(double *tmp, double *tmpV, int n){
    for (int i=0; i<n; i++)
        if (tmp[i] > angleMax) {  tmp[i] = angleMax;  tmpV[i] = -tmpV[i]; }
        else if (tmp[i] < angleMin)   {    tmp[i] = angleMin;    tmpV[i] = -tmpV[i]; }
    //  if the particle reach the boudary, v need to inverse
}

void checkOrigin(double *tmp, double *tmpV, int n){
    for (int i=0; i<n; i++)
        if (tmp[i] > angleMax) {  tmp[i] = angleMax;  tmpV[i] = -tmpV[i]; }
        else if (tmp[i] < angleMin)   {    tmp[i] = angleMin;    tmpV[i] = -tmpV[i]; }
    //  if the particle reach the boudary, v need to inverse
}

double getInertiaWeight(double it, double MaxIt){
    return wMax- (wMax-wMin)*(it+1)/MaxIt;
}

int rouletteWheel(myRep & rep){
    int sz = static_cast<int>(rep.size());
    return rand() % sz;
}


// SRD: sqrt(rep[i].cost - particle.cost), minSRD is the best rep member for particle
int getGBest(Particle &particle, myRep &rep) {
    int size = static_cast<int>(rep.size());

    double **f;            //  前size个是rep中元素的cost，最后一个放当前particle的cost  f[size+1][objectiveNumber]
    f = new double *[size+1];
    for (int i = 0; i < size+1; i++)
        f[i] = new double [objectiveNumber];

    myRep ::iterator a = rep.begin();
    for (int i = 0; i < size; i++){
        for (int j = 0; j < objectiveNumber; j++)
            f[i][j] = a->Cost[j];
        a++;
    }
    for (int j = 0; j < objectiveNumber; j++)
        f[size][j] = particle.Cost[j];

    unitFunction(f, size+1, objectiveNumber);

    double minSRD = 1000;
    int ans = 0;

    for (int i = 0; i < size; i++){
        double tmp = getSRD(f, i, size);
        if ( tmp < minSRD){
            ans = i;
            minSRD = tmp;
        }
    }

    freeUpSpace(f, size+1);
    return  ans;
}

void unitFunction(double **f, int n, int m){
    double maxCost[m], minCost[m];      // minCost[0] : the min value of the first energy function
    for (int i = 0; i < m; i++){
        maxCost[i] = -9999999999;
        minCost[i] = std::numeric_limits<double>::max();
    }

    for (int i = 0; i < n; i++)
        for (int j = 0; j< m; j++) {
            if (f[i][j] > maxCost[j])
                maxCost[j] = f[i][j];
            if (f[i][j] < minCost[j])
                minCost[j] = f[i][j];
        }

    /*
     * calculate newCost for every f
     */
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            f[i][j] = (f[i][j] - minCost[j]) / (maxCost[j] - minCost[j]);
}

bool comp_deleteArray(const int &a, const int &b){
    return a > b;
}

double getSRD(double **f, int k, int size) {
    double ans = 0;
    for (int i = 0; i < objectiveNumber-1; i++){
        ans += sqrt(fabs(f[size][i] - f[k][i]));        //  f[size][i] is the particle's cost
    }
    return ans;
}

double dis1(double x, double y, double z){
    double ans = (x-1)*(x-1) + (y-5)*(y-5) + (z+12)*(z+12);
    return sqrt(ans);
}

double dis2(double x, double y, double z){
    double ans = (x+1)*(x+1) + (y-3)*(y-3) + (z+7)*(z+7);
    return sqrt(ans);
}

char * getRelativeAA(char abbreviation){
    char *str;
    str = new char [3];
    switch (abbreviation){
        case 'A':   strcpy(str, "ALA"); break;
        case 'R':   strcpy(str, "ARG"); break;
        case 'N':   strcpy(str, "ASN"); break;
        case 'D':   strcpy(str, "ASP"); break;
        case 'C':   strcpy(str, "CYS"); break;
        case 'Q':   strcpy(str, "GLN"); break;
        case 'E':   strcpy(str, "GLU"); break;
        case 'G':   strcpy(str, "GLY"); break;
        case 'H':   strcpy(str, "HIS"); break;
        case 'I':   strcpy(str, "ILE"); break;
        case 'L':   strcpy(str, "LEU"); break;
        case 'K':   strcpy(str, "LYS"); break;
        case 'M':   strcpy(str, "MET"); break;
        case 'F':   strcpy(str, "PHE"); break;
        case 'P':   strcpy(str, "PRO"); break;
        case 'S':   strcpy(str, "SER"); break;
        case 'T':   strcpy(str, "THR"); break;
        case 'W':   strcpy(str, "TRP"); break;
        case 'Y':   strcpy(str, "TYR"); break;
        default :   strcpy(str, "VAL"); break;
    }
    return str;
}



void strReplace(char *src, const char *shortStr, const int start, const int length){
    char *newSrc = new char [256]; //strlen(src) - length + strlen(shortStr)
    strncpy(newSrc, src, start);
    newSrc[start] = '\0';
    cout << "newSrc=" << newSrc << endl;
    strcat(newSrc, shortStr);
    strcat(newSrc, src+start+length);
    strcpy(src, newSrc);
    delete newSrc;
}

int strToInt(char *str){
    int ans = 0;
    int len = strlen(str);
    for (int i=0; i<len; i++){
        ans = ans*10 + str[i] - '0';
    }
    return ans;
}

/*
 [argv] <0>       <1>                           <2>                     <3>                 <4>                     <5.6>
     ./AIR /home/advanced/Data/TR829    7                            3000                 50                    1 2
       <input address>                 <standard particle numbers>   <iteration times>    <total population>    <2 reps for TM_align>
*/
int getInputParameter(char **argv) {
    inputSize = strToInt(argv[2]);
    MaxIt = strToInt(argv[3]);
    Population = strToInt(argv[4]);
}

int getMax(int a, int b){
    if ( a > b )return a;
    else return b;
}

bool doubleEqual(double a, double b){
    if (abs(a-b) < Criterion)
        return 1;
    else
        return 0;
}

//  help to print parameters
void printPara(char *s, int k){
    printf("%s : %d\n", s, k);
}

void printPara(char *s, double d){
    cout << s << " : " << d << endl;
}
void printPara(char *s, char *str){
    printf("%s : %s\n", s, str);
}
