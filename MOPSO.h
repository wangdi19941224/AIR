#ifndef MOPSO_H
#define MOPSO_H

#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdio>
#include <cstring>
#include <list>
#include <cstdlib>
#include <ctime>
#include <pthread.h>
#include <vector>
#include <iomanip>
#include <stdexcept>
#include <limits>
#include <algorithm>
#include <string>

using namespace std;

struct structBest;
struct Particle;
struct Rep;
struct Atom;
struct POINT;
struct VECTOR;
struct ANGLE;
struct EnergyValueStruct;
struct REP_COUNT;
template <typename _MD>
struct commonStruct;
template <typename _MD>
class Matrix;


typedef list<Rep> myRep;            //   a wonderful list, the data struct of the repositpry, which
typedef const char CONCH;
//store the best solutions that have been found (pateto front)


// fileDisposal
    // Input Data

int getLines(char *address);
void getAtom(Atom *atom, int cntLines, int &atomNum, char *address);        // find the word "ATOM" in the file

char * inputSeq();
void inputParticle(Particle &particle, int k, char *seq);
void inputParticles(Particle *particle);
void inputParOrigin(Particle &particle, int fileNum, int cntLines);
void inputParPhi(Particle &particle, int fileNum);
    // a function that will open file and throw error if it can't open file
void openFile(const char fileDir[10000], ifstream &infile);
    // a function that will remove file and throw error if it can't remove file
void removeFile(const char * str);
void removeMultiplyFile_afterInput();

    // for the object function
void printPdb(Particle * particle);
void printPdb(Particle * particle, const char *fileName);
void printCopiedOrigin(const Particle &particle, const int repNum);
void printParticle(Particle * particle);
void printPdb(list<Rep>::iterator it, const int &repNum);
    // for debug
        // print rep in several files
void printRep(myRep & rep, int len);
        // print the particles' status in one file
void printParticleCost(Particle * particle, int iterator);
    // print answer rep
    void outputAnswer(myRep &rep);
    //  print run time
    void printTime(const int choice, const int loopTimes);

// print time for ever different energe function
void printCostTime(long rosettaTime, long calTime , long charmmTime, long mybinTime);
    //  create a new fold for the answer
void createNewFold();

//  multiply better input
void multiplyBetterInput();
void multiplyParticle(Particle &particle, int repNum);

void runScwrl(int index);                   //  run scwrl program after print .pdb
void removeTmp(int index);

// get the cost for the object function 1 and 2
double getCost_rosetta(int i);
double getCost_calRWplus(int i);
double getCost_charmm(int i);
double getCost_mybin(int i);
void rewriteTemplate(int i, const char *templateFile, const char *tmpIFile, const char *tmpI);
double getEnergyValue(const char *fileName);
void printCurrentId(int id);

// initialize particle
void initializeParticles(Particle *particle);
void initializeParticle(Particle &particle);

// put input particles into rep
void updateRep(Particle *particle, myRep &rep, int loopTimes);

//  MOPSOFunction
    //  apply the formulation of PSO
void PSOAdaptionForPhi(Particle &particle, myRep &rep, int it, const int i);

    //  decide if particle[i] is dominated, for all i
void decideDominated(Particle * particle);

    //  compare the new position with the pBest and update it
void updatePBest(Particle * particle);

    // calculate f(x)
void getAllParticleCostForTest(Particle * particle);
void getAllParticleCost(Particle * particle);

void *getAParticleCost(void* particle);

    //if the result is true, then cost1 is dominated by cost2
bool isDominated(double *cost1, double *cost2);

    //    add the new non-dominated particles into rep / eliminate some in the rep that is not so good now
void putNewParticleIntoRep(Particle * particle, myRep & rep, int iterator);

    //decide whether particle[i] can be put into rep or not, as well as eliminate the rep particles
    //which are dominated by particle[i]
    // if return true then can put
bool canPutIntoRep(double *Cost, myRep & rep);

void sieveRep(Particle *particle, myRep &rep, const int &loopTimes);         //  if there are too many paticles in rep, restrict it in a set value


//  MOPSOAidFunction
    //  cat c1-ss-c2  e.g. :  "mj" , 520 , "forever"  ->  mj520forever        used to conbine file name
char *catStrIntStr(const char *c1, int ss, const char *c2);
char *catStrIntStr(const char *c1, const char *c2, int ss, const char *c3);

    //  cat c1-c2
char *catStrStr(CONCH *c1, CONCH *c2);
char *catStrStr(CONCH *c1, CONCH *c2, CONCH *c3);
char *catStrStr(CONCH *c1, CONCH *c2, CONCH *c3, CONCH *c4);
char *catStrStr(CONCH *c1, CONCH *c2, CONCH *c3, CONCH *c4, CONCH *c5);
char *catStrStr(CONCH *c1, CONCH *c2, CONCH *c3, CONCH *c4, CONCH *c5, CONCH *c6);
char *catStrStr(CONCH *c1, CONCH *c2, CONCH *c3, CONCH *c4, CONCH *c5, CONCH *c6, CONCH *c7);
char *catStrStr(CONCH *c1, CONCH *c2, CONCH *c3, CONCH *c4, CONCH *c5, CONCH *c6, CONCH *c7, CONCH *c8);

    //  copy p2 to p1
void cpyDoubleArray(double *& p1, double *& p2, int n);
void cpyIntArray(int *& p1, int *& p2, int n);

    //  check if V is out of range
void checkV(double *tmp, int n);

//  check if position is out of range
void checkOrigin(double *tmp, double *tmpV, int n);

    //  check if position is out of range
void checkP(double *tmp, double *tmpV, int n);

    //  get the right rep_h for this update time
int rouletteWheel(myRep & rep);

// get GBest through SRD
int getGBest(Particle &particle, myRep &rep);

// get SRD
double getSRD(double **f, int k, int size);

// unitize 归一化
void unitFunction(double **f, int n, int m);

bool comp_deleteArray(const int &a, const int &b);

    //  get the w
double getInertiaWeight(double it,double MaxIt);

    //  for test object function
double dis1(double x, double y, double z);
double dis2(double x, double y, double z);

    //  get the AA that pr represent
char * getRelativeAA(char pr);

void strReplace(char *src, const char *shortStr, const int start, const int length);

int strToInt(char *str);

int getInputParameter(char **argv);

int getMax(int a, int b);

bool doubleEqual(double a, double b);

//  help to print parameters
void printPara(char *s, int k);
void printPara(char *s, double d);
void printPara(char *s, char *str);

//  convert rotation to coordinary
void convertRotationToCoordinary(Particle &particle);


//  Matrix Function
    //  ans = a*b
void matrixProduct(double a[4][4], double b[4][4], double ans[4][4]);

    //  set a[4][4] to be full of 0
void initialize_0(double a[4][4]);

    //  dan wei(chinese) matrix
void initialize_1(double a[4][4]);

    //  get the R in the paper
void getR(int it, double * a, double theta, double R[4][4]);

    //  get the new coordinate for itAtom
void getNewCoordinary(int itAtom, double R[4][4], double * a);

void cpyOriginTOAddO_origin(int it, double * x, double * xAddO);

void getNewCoordinaryForO(int it, double R[4][4], double * xAddO, int numAtom);

    //  newR = preR
void cpyMatrix(double preR[4][4], double newR[4][4]);

    //   dan wei hua (chinese)
void unitization(double v[3]);

    //set Q[3][3] the paper used
void getQ(double Q[4][4], double v[3], double theta);


//  get the best reps
void sortRepByAngle(myRep &rep, ANGLE *angle);
void sortRepByLambda(myRep &rep);
void getCostInUnit(myRep &rep, double **costInUnit);
void getAngle(int i, const POINT *point, ANGLE *angle);
double vectorDot(VECTOR u, VECTOR v);
double cross(VECTOR u, VECTOR v);
double length(VECTOR u);
bool pdPlusPi(POINT first, POINT second, POINT third);
void angleQsort(int l, int r, ANGLE *a);
void energyQsort(int l, int r, EnergyValueStruct *a);
void pointQsort(int l, int r, POINT *a);

// setTime
void setTime();
void getEndTime();

// dispose the input parameters and files previously
void preDisposeInputParametersAndFiles(char **argv);
void disposePDB();
void createSeqTxt();
void createParticleTxt(int k);
void createPhi(int k);
bool isImportantAtom(const Atom atom);
char getAbbreviation(char *str);
void getArgv();
void getParameters_for_TMalign();

// check similarity
void checkAllParticleSimilarity(Particle *particle);
void *checkOneParticleSimilarity(void *p);
void checkOneParticleSimilarity(Particle &particle);
void runTM_score(int index);                //  run TM_score to get the similarity of 2 particle
double getTM_score(int index);
bool particleSimilarity(double TM_scrore);      // if return 1 then two particles are similar
void becomeInitialParticle(Particle &particle);

//free up space
void freeUpSpace(Particle *& particle);
void freeUpSpace(double **f, int rowsize);

// apply variables
void applyVariable();
// release space
void releaseSpace(double **x, const int &n);
void releaseSpace(double *x);

// code for Debug
void printIn(int option);
void printOut(int option);

// disposeTMScore
/*
 * Position 0 and Position 1 refer to first AA while Position 1 is always 180
 * Position 2 ,Position 3, Position 4 refer to second[2] AA while Position 4 is always 180
 * Position 5 ,Position 6, Position 7 refer to third[3] AA while Position 7 is always 180
 * Position 8 ,Position 9, Position 10 refer to fourth[4] AA while Position 10 is always 180
 * ...  if there's only 5 AA, then Position 11 refer to the last AA
 * so for not first nor last AA, (Kth_AA - 2) * 3 + 2, (Kth_AA - 2) * 3 + 3, (Kth_AA - 2) * 3 + 4 is Mth_Atom
 */
void updateVelocityCheck(const int &iterationTimes);
void setVelMax(double *VelMax, const int &seq_AANum, const int &numAA, const int &VelMax_Len);
int getVelMax_by_TMalign(const int &firRep_for_TMalign, const int &secRep_for_TMalign);

/*
 * Matrix function
 */
/*
     * Sum of two Matrix
     */
template <typename _MD>
Matrix<_MD> operator+(const Matrix<_MD> &a, const Matrix<_MD> &b);

template <typename _MD>
Matrix<_MD> operator+(const Matrix<_MD> &a, const _MD &b);

template <typename _MD>
Matrix<_MD> operator-(const Matrix<_MD> &a, const Matrix<_MD> &b);

Matrix<double> operator-(const Matrix<double> &a, const Matrix<int> &b);

Matrix<double> operator-(const Matrix<int> &a, const Matrix<double> &b);

Matrix<double> operator-(const Matrix<int> &y, const double d);

Matrix<double> operator-(const Matrix<double> &y, const int d);

Matrix<double> operator-(const double d, const Matrix<int> &y);

Matrix<double> operator-(const int d, const Matrix<double> &y);

template <typename _MD>
Matrix<_MD> operator-(const Matrix<_MD> &y, const _MD d);

template <typename _MD>
Matrix<_MD> operator-(const _MD d, const Matrix<_MD> &y);

template <typename _MD>
bool operator==(const Matrix<_MD> &a, const Matrix<_MD> &b);

template <typename _MD>
Matrix<_MD> operator-(const Matrix<_MD> &a);

/*template <typename _MD>
Matrix<_MD> operator-(Matrix<_MD> &&a);*/
/*
 * multiplication of two matrix
 */
template <typename _MD>
Matrix<_MD> operator*(const Matrix<_MD> &a, const Matrix<_MD> &b);
/*
 * matrix multiply a number
 */
template <typename _MD>
Matrix<_MD> operator*(const Matrix<_MD> &a, const _MD &b);

template <typename _MD>
Matrix<_MD> operator*(const _MD &b, const Matrix<_MD> &a);

template <typename _MD>
Matrix<_MD> operator/(const Matrix<_MD> &a, const _MD &b);

Matrix<double> operator/(const Matrix<double> &a, const int &b);

Matrix<double> operator/(const Matrix<int> &a, const double &b);

/*
 * if a[i][j] > b[i][j], c[i][j] = 1, or c[i][j] = 0
 */
template <typename _MD>
Matrix<int> operator>(const Matrix<_MD> &a, const Matrix<_MD> &b);


/*
 * if a[i][j] <= b[i][j], c[i][j] = 1, or c[i][j] = 0
 */
template <typename _MD>
Matrix<int> operator<=(const Matrix<_MD> &a, const Matrix<_MD> &b);

/*   .*   规模相同的矩阵对应位置数相乘    mulCorNum(M,M)  //  multiply corresponding number
     ./  				             divCorNum(M,M)  //  two matrix with the same size
     .^				                 powCorNum(M,D)  */
template <typename _MD>
Matrix<_MD> mulCorNum(const Matrix<_MD> &a, const Matrix<_MD> &b);

Matrix<double> mulCorNum(const Matrix<double> &a, const Matrix<int> &b);

Matrix<double> mulCorNum(const Matrix<int> &a, const Matrix<double> &b);

template <typename _MD>
Matrix<_MD> divCorNum(const Matrix<_MD> &a, const Matrix<_MD> &b);

Matrix<double> divCorNum(const Matrix<double> &a, const Matrix<int> &b);

Matrix<double> divCorNum(const Matrix<int> &a, const Matrix<double> &b);

template <typename _MD>
Matrix<_MD> powCorNum(const Matrix<_MD> &a, const int &p);

/*
 * y = [1,2,3;4,5,6;7,8,9];
 * y(1,1:3) = [5,6];
 */
template <typename _MD>
Matrix<_MD> getArrayFromMatrix(const Matrix<_MD> &y, const int refRow, const size_t fromCol, const size_t endCol);

Matrix<int> initArray(const int &startNum, const int &delta, const int &endNum);

template <typename _MD>         // repmat(MATRIX, n, m) /  repmat(double, n, m)
Matrix<_MD> repmat(const Matrix<_MD> &y, int n, int m);

// repmat(MATRIX, n, m) /  repmat(double, n, m)
Matrix<double> repmat(const double d, int n, int m);

Matrix<int> repmat(const int d, int n, int m);

template <typename _MD>
Matrix<_MD> transpose(const Matrix<_MD> &a);

//  unit MATRIX
template <typename _MD>
Matrix<_MD> I(const size_t &n);

template <typename _MD>
Matrix<_MD> matPow(Matrix<_MD> mat, int b);

/*
 * Matrix mat;
 * mat^(1/2)
 */
template <typename _MD>
Matrix<_MD> matRoot(const Matrix<_MD> &mat, double p);

/*  matrixSort : sort by column
 *  [ig,p] = sort(rp);    p : the same size with rp, p[i][j] : p[i][j]是rp矩阵j列的第几个
 */
template <typename _MD>
Matrix<int> matrixColSort(const Matrix<_MD> y);

template <typename _MD>
void colSort(int l, int r, commonStruct<_MD> *a);

/*
 * MATRIX x, MATRIX y , z = x(y)
 * the same size with y, z[i][j] is the y[i][j](th) number in x Matrix(column first)
 */
template <typename _MD>
Matrix<_MD> matBracket(const Matrix<_MD> &x, const Matrix<int> &y);

/*
 * t = [1,2,3; 4,5,6];  j = [2,2,2];
 * t(0,j) = [3,3,3];
 * t(1,j) = [6,6,6]
 * this function is to obtain t(i,j); when j is a matrix
 */
template <typename _MD>
Matrix<_MD> matRowBracket(const Matrix<_MD> &x, const int &i, const Matrix<int> &y);

template <typename _MD>
_MD getKthNumber(const Matrix<_MD> &mat, int k);

Matrix<double> randMatrix(const int n, const int m);

/*
 * change 1 to 0 / 0 to 1
 */
template <typename _MD>
Matrix<_MD> NOT(const Matrix<_MD> &y);

template <typename _MD>
std::ostream &operator<<(std::ostream &os, const Matrix<_MD> &mat);

template <typename _MD>
Matrix<_MD> inputMatrix(const size_t &n, const size_t &m, const _MD &type);

/*
 * randFixedSum
 */
template <typename _XD>
Matrix<double> randFixedSum(int n, int m, _XD S, _XD a, _XD b);

//  const statement
extern const char rootAddress[];
extern char *inputAddress;
extern char *logAddress;
extern char *energyOutputAddress;
extern double *VelMax;
extern const char *energyFileAddress;
extern const char *tempFileAddress;
extern const char *mybinAddress;
extern const char *defaultFileAddress;
extern const char *QUACKoutFileAddress;
extern const char *calFileAddress;
extern const char *charmmFileAddress;
extern const char *refine_1Address;
extern const char *TM_scoreAddress;
extern char *answerAddress;
extern const char *draftAddress;
extern const char *scoreAddress;
extern const char *databaseAddress;
extern const char *calRWplusAddress;
extern const char *strideAddress;
extern const char *TM_alignAddress;

extern const double angleMin;
extern const double angleMax ;


// MOPSO Settings
extern int inputSize;                  // Population Size
extern const int nRep;                // Repository Size
extern int Population;
extern int  MaxIt;          // Maximum Number of Iterations
extern int multiplyNumber; // copy the best input particle
extern const double Criterion;
extern const int lambdaLoopTimes;
extern const int tidSize;
extern time_t starTime;
extern time_t endTime;
extern char **argv;
extern int startNum;
extern int numAA;
extern int currentIteration_Times;
extern int times_for_each_play;
extern double VEL_SMALL_RANGE;
extern double VEL_BIG_RANGE;

extern const double phi1;
extern const double phi2;
extern const double phi;
extern const double chi ;    // 0.73
extern const int bufferLen;

extern const double wMin;                        //  Inertia Weight
extern const double wMax;
extern const double wDamp;                       //  Inertia Weight Damping Ratio

//extern const double c1;                 //  Personal Learning Coefficient
//extern const double c2;                 //  Global Learning Coefficient
extern const double c1max;
extern const double c1min;
extern const double c2max;
extern const double c2min;

extern const int objectiveNumber;      //  Multiple Objectives
extern const double TM_scoreThreshold; // the threshold of TM-score
extern const double PI;
extern const double INF;
extern const double realMax;
extern const double tiny;

// my tools
extern Particle *particle;
extern int *par_firNumber;
extern int *par_secNumber;
extern myRep rep;
extern ANGLE *angle;
extern int *sortAns;
extern REP_COUNT *rep_count;

//  struct definition
// pBest
struct structBest{
    double *Position;
    double *addO_origin;
    double *Cost;
};

struct Particle{
    double * origin;                          //    the coordinate(x,y,z) of every atom (not include atom O),
    //  three indexes represent an Atom's position
    double * addO_origin;                   // the coordinate(x,y,z) of every atom (include the atom O)
    double * Position;                      //  the angles
    double * old_position;                  //  used in the calculation of angle change
    double * Velocity;                      //  will be used in the MOPSO
    double * Cost ;                         //  will be used int the MOPSO, store every f(x) the particle have
    int * GridIndex, * GridSubIndex;
    bool dominated;                     //   represent if the particle is dominated by some other,
    //  which means it is certainly not the best answer
    int numAA, numAtom;         //  the number of AA and the number of atoms
    int sizeOfOrigin, sizeOfAddO_origin, sizeOfPosition, sizeOfVelocity;    //  the size of arrays
    int index;
    char *seq;
    int startNubmber;
    structBest Best;                    //  pBest
    list<structBest> BestArray;
};

struct Rep{
    double *Cost;
    double *addO_origin;
    double *Position;
    int sizeOfAddO_origin;
    int iterator;
    char *seq;
    int startNumber;
    int numAA;
};

struct Atom {
    char name[10], group[10];
    int firNumber, secNumber;
    double x, y, z;
};

struct POINT{
    double x, y;
    int id;
};

struct VECTOR{
    double x, y;
    VECTOR(double x1, double y1){
        x = x1;
        y = y1;
    }
};

struct ANGLE{
    double value;
    int id;
};

struct EnergyValueStruct{
    double value;
    int id;
};

template <typename _MD>
struct commonStruct{
    _MD value;
    int id;
};

struct REP_COUNT{
    int id, amount;
    bool operator<(const REP_COUNT &rhs) { return amount > rhs.amount;}
};

template <typename _MD>
class Matrix{
protected:
    size_t n_rows;
    size_t n_cols;
    std::vector<std::vector<_MD> > data;
    class RowProxy {
    private:
        std::vector<_MD> &row;                                   //***** & quote, give a variable another name
    public:
        RowProxy(std::vector<_MD> &_row) : row(_row) {}         // ***** &
        _MD &operator[](const size_t &pos){
            return row[pos];
        }
    };
    class ConstRowProxy {
    private:
        const std::vector<_MD> &row;
    public:
        ConstRowProxy(const std::vector<_MD> &_row) : row(_row) {}
        const _MD &operator[](const size_t &pos) const {
            return row[pos];
        }
    };

public:
    Matrix() {};
    Matrix(const size_t &_n_rows, const size_t &_n_cols)
            : n_rows(_n_rows), n_cols(_n_cols), data(std::vector<std::vector<_MD> >(_n_rows, std::vector<_MD>(n_cols))) {}
    Matrix(const size_t &_n_rows, const size_t &_n_cols, _MD fillValue)
            : n_rows(_n_rows), n_cols(_n_cols), data(std::vector<std::vector<_MD> >(_n_rows, std::vector<_MD>(_n_cols, fillValue))) {}
    Matrix(const Matrix<_MD> &mat)
            : n_rows(mat.n_rows), n_cols(mat.n_cols), data(mat.data) {}
//    Matrix(Matrix<_MD> &&mat) noexcept
//            : n_rows(mat.n_rows), n_cols(mat.n_cols), data(mat.data) {}
    Matrix<_MD> & operator=(const Matrix<_MD> &mat){
        n_rows = mat.n_rows;
        n_cols = mat.n_cols;
        data = mat.data;
        return *this;
    }
//    Matrix<_MD> &operator=(Matrix<_MD> &&mat){
//        n_rows = mat.n_rows;
//        n_cols = mat.n_cols;
//        data = mat.data;
//        return *this;
//    }
    inline const size_t & RowSize() const {
        return n_rows;
    }
    inline const size_t & ColSize() const {
        return n_cols;
    }
    RowProxy operator[](const size_t &Kth){
        return RowProxy(data[Kth]);
    }

    const ConstRowProxy operator[](const size_t &Kth) const {
        return ConstRowProxy(data[Kth]);
    }

    /*
 * w(i,2:i+1) = tmp1 + tmp2;
 */
    Matrix<_MD> matRowChange(const int &row, const int &startCol,
                             const int &endCol, const Matrix<_MD> &x);
    ~Matrix() {};
};



#endif // MOPSO_H
