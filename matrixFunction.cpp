#include "MOPSO.h"

void matrixProduct(double a[4][4], double b[4][4], double ans[4][4]){
    initialize_0(ans);
    for (int i=0; i<4; i++)
        for (int j=0; j<4; j++)
            for (int k=0; k<4; k++)
                ans[i][j] += a[i][k] * b[k][j];
}

void initialize_0(double a[4][4]){
    for (int i=0; i<4; i++)
        for (int j=0; j<4; j++)
            a[i][j] = 0;
}

void initialize_1(double a[4][4]){
    for (int i=0; i<4; i++)
        for (int j=0; j<4; j++)
            if (i==j)   a[i][j] = 1;
            else a[i][j] = 0;
}

//  a is the coordinary array,  theta is the rotation angle
void getR(int itAtom, double * a, double theta, double R[4][4]){
    double v[3];
    int ii=(itAtom-2) *3, jj=(itAtom-1) *3;       // i: the x index of (itAtom-2)
    for (int k=0; k<3; k++)
        v[k] = a[ii+k] - a[jj+k];                             // v = Q(i-1) - Qi
    unitization(v);                                         // v[0]^2 + v[1]^2 + v[2]^2 = 1
    double Q[4][4];
    initialize_0(Q);
    getQ(Q, v, theta);          //set Q[3][3]
    for (int i=0; i<4; i++)
        for (int j=0; j<4; j++)
            if (i==3 || j==3)
                Q[i][j] = 0;
    Q[3][3] = 1;

    double T[4][4], transposeT[4][4];
    initialize_1(T);
    initialize_1(transposeT);
    jj=(itAtom-1) *3;
    for (int p=0; p<3; p++){
        T[p][3] = a[jj+p];
        transposeT[p][3] = -a[jj+p];
    }
    double tmp[4][4];
    matrixProduct(T, Q, tmp);
    matrixProduct(tmp, transposeT, R);
}

void getQ(double Q[4][4], double v[3], double theta){
    theta = theta/180  * M_PI;
    double c = cos(theta), s = sin(theta);
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            switch ( (j+3-i) % 3) {
                case 0 :
                    Q[i][j] = v[i]*v[i]+ (1-v[i]*v[i])*c;
                    break;
                case 1 :
                    Q[i][j] = v[(i+1)%3] * v[i] * (1-c) - v[(i+2) % 3]*s;
                    break;
                case 2 :
                    Q[i][j] = v[(i+2) %3] * v[i] *(1-c) + v[(i+1) % 3]*s;
                    break;
                default: break;
            }
        }
    }
}

void unitization(double v[3]){
    double ans=0;
    for (int i=0; i<3; i++)
        ans += v[i] * v[i];
    ans = sqrt(ans);
    for (int i=0; i<3; i++)
        v[i] /= ans;
}

void getNewCoordinary(int itAtom, double R[4][4], double * a){
    double b[4];
    for (int i=0; i<3; i++)
        b[i] = a[ itAtom*3+i ];
    b[3] = 1;
    for (int i=itAtom*3; i<itAtom*3+3; i++){
        a[i] = 0;
        for (int k=0; k<4; k++)
            a[i] += R[i-itAtom*3][k] * b[k];
    }
}

void cpyOriginTOAddO_origin(int itAtom, double * x, double * xAddO){
    int i1 = itAtom*3;                          // the begin position of this AA
    int i2;
    if (itAtom==3)  i2 =4*3;                    // the end sequence of this AA
    else {
        int numOfAAPre = itAtom/3;
        int numOfAtom = numOfAAPre*4 + itAtom % 3;
        i2 = numOfAtom * 3;
    }
    for (int k=0; k<3; k++)
        xAddO[i2+k] = x[i1+k];
}

void getNewCoordinaryForO(int itAtom, double R[4][4], double * xAddO, int numAtom){
    if (itAtom==numAtom-1) itAtom += 1;
    int numOfAA = itAtom / 3;
    int numOfAtomAddO = numOfAA*4-1;
    getNewCoordinary(numOfAtomAddO, R, xAddO);
}

void cpyMatrix(double preR[4][4], double newR[4][4]){
    for (int i=0; i<4; i++)
        for (int j=0; j<4; j++)
            preR[i][j] = newR[i][j];
}

