//
// Created by Zyj on 2018/1/21.
//

#include "DH_rotation.h"
COORDINATES **frame;
COORDINATES **frame_ForO;
COORDINATES **frame_ForC;

/*  create frame for each particles' each atom
  set DH parameter
  init frame.
  It's very difficult to calculate coordinates of O, because O is not included in the chain
  so we need to do some extra work and write some new functions to get the DH parameter of O
  */
void initDH_rotation(){
    frame = new COORDINATES *[Population];

    for (int i = 0; i < Population; i++) {
        int numAtom = particle[i].numAtom;
        frame[i] = new COORDINATES[numAtom + 2];        //  add two in order to give world frame and help calculate the last atom's DHparameter
        /*
         * frame[i][0] refers to world coordinate(0, 0, 0)/ (1 0 0) (0 1 0) (0 0 1)
         * so frame[i][1] refers to the first atom in the particle
         */
        frame[i][0].origin = POINT3d(0, 0, 0);

        for (int j = 0; j < numAtom; j++){
            frame[i][j+1].origin.x = particle[i].origin[j * 3 + 0];
            frame[i][j+1].origin.y = particle[i].origin[j * 3 + 1];
            frame[i][j+1].origin.z = particle[i].origin[j * 3 + 2];
        }

        frame[i][numAtom + 1].origin = POINT3d(0, 0, 0);            // in order to help establish last atom's frame
    }
    createChainCoordinates();
    setDHparameter();

    initDH_paraForO();
}

void convertRotationToOrigin(Particle &particle){
    double transferMat[4][4];
    initialize_1(transferMat);
    int numAtom = particle.numAtom;
    int pID = particle.index;

    /*
     * first rotation angle is not included in phi.txt
     * so the first two atom : N & C 's origin can't change
     * because the key rotates along the key's direction, so the second C's origin also not change
     * origin change will begin from the fourth Atom N, whose atomID is 3
     * we need to prepare transfer matrix for it now
     * frame[0] refers to the world coordinate in DH method, so frame[1] refers to atomID = 0
     * numAtom have ignored atom O, from the second N
     */

    for (int atomID = 0; atomID < numAtom; atomID++){
        int frameID = atomID + 1;
        double theta = getTheta(particle, atomID);

        frame[pID][frameID].theta = theta;
        updateDHmatrix(frame[pID][frameID]);
        double tmpAns[4][4];
        matrixProduct(transferMat, frame[pID][frameID].DH_matrix, tmpAns);
        copyMatrix(tmpAns, transferMat);
        updateParticleCoordinate(particle, atomID, transferMat);
        cpyOriginTOAddO_origin(atomID, particle.origin, particle.addO_origin);

        // update atom O, when we meet CA
        if (atomID % 3 == 1){
            getNewCoordinateForO(particle, transferMat, atomID);
        }
    }
}

/*
 * first we need to calculate transfer matrix of C, then matrix of O
 */
void getNewCoordinateForO(Particle &particle, double transferMat[4][4], const int &atomID) {
    int AA_Id = atomID / 3;
    double tmpAns[4][4], tmpAns_C[4][4];
    int pID = particle.index;

    // update C's transfer matrix, current atomID is for CA
    double theta_ForC = getTheta(particle, atomID+1);
    frame_ForC[pID][AA_Id].theta = theta_ForC;
    updateDHmatrix(frame_ForC[pID][AA_Id]);
    matrixProduct(transferMat, frame_ForC[pID][AA_Id].DH_matrix, tmpAns_C);

    //  update O's transfer matrix and calculate its new coordinary
    matrixProduct(tmpAns_C, frame_ForO[pID][AA_Id].DH_matrix, tmpAns);
    // calculate O's new coordinate
    VECTOR3d B(0, 0, 0);
    VECTOR3d ans = tmpAns * B;

    //printf("ans : %lf  %lf  %lf\n\n", ans.x, ans.y, ans.z);

    particle.addO_origin[AA_Id * 12 + 9 ] = ans.x;
    particle.addO_origin[AA_Id * 12 + 10] = ans.y;
    particle.addO_origin[AA_Id * 12 + 11] = ans.z;
}

void createChainCoordinates(){
    for (int i = 0; i < Population; i++){
        int numAtom = particle[i].numAtom;

        frame[i][0].X = VECTOR3d(1, 0, 0);
        frame[i][0].Y = VECTOR3d(0, 1, 0);
        frame[i][0].Z = VECTOR3d(0, 0, 1);
        /*
         * frame[i][0] refers to world coordinate(0, 0, 0)/ (1 0 0) (0 1 0) (0 0 1)
         * so frame[i][1] refers to the first atom in the particle
         * atom in particle refers to index [1]~[numAtom] in frame
         */
        for (int j = 1; j <= numAtom; j++){
            frame[i][j].Z = std_VECTOR3d(frame[i][j + 1].origin - frame[i][j].origin) ;
            frame[i][j].X = std_VECTOR3d(cross_VECTIR3d(frame[i][j].Z, frame[i][j-1].Z));
            frame[i][j].Y = std_VECTOR3d(cross_VECTIR3d(frame[i][j].Z, frame[i][j].X));
        }
    }

}

void setDHparameter(){
    for (int i = 0; i < Population; i++){
        int numAtom = particle[i].numAtom;

        for (int j = 1; j <= numAtom; j++){

            frame[i][j].d = dist_POINT3d(frame[i][j-1].origin, frame[i][j].origin);
            frame[i][j].a = 0;
            frame[i][j].init_theta = findAngle(frame[i][j-1].X, frame[i][j].X, frame[i][j-1].Z);
            frame[i][j].alpha      = findAngle(frame[i][j-1].Z, frame[i][j].Z, frame[i][j].X);
        }
    }

}
/*
 * It's very difficult to calculate coordinates of O, because O is not included in the chain
 * so we need to do some extra work and write some new functions to get the DH parameter of O
 * I establish 2 new array : frame_ForO, frame_ForC
 * if we can define frame_ForO at the beginning, then we can calculate new origin of O in each updation, from C to O
 * to define frame_ForO, we need to know the coordinate of atom before O, which is C
 *
*/
void initDH_paraForO(){
    frame_ForO = new COORDINATES *[Population];
    frame_ForC = new COORDINATES *[Population];
    // init coordinate
    for (int i = 0; i < Population; i++){
        int numAA = particle[i].numAA;                  // one AA has one O
        frame_ForO[i] = new COORDINATES[numAA];
        frame_ForC[i] = new COORDINATES[numAA];

        // get C and O's coordinate in an AA
        for (int j = 0; j < numAA; j++){
            frame_ForC[i][j].origin.x = frame[i][j * 3 + 3].origin.x;
            frame_ForC[i][j].origin.y = frame[i][j * 3 + 3].origin.y;
            frame_ForC[i][j].origin.z = frame[i][j * 3 + 3].origin.z;

            frame_ForO[i][j].origin.x = particle[i].addO_origin[j * 12 + 9];
            frame_ForO[i][j].origin.y = particle[i].addO_origin[j * 12 + 10];
            frame_ForO[i][j].origin.z = particle[i].addO_origin[j * 12 + 11];
        }
    }
    createFrameXYZ_ForO();
    setDH_parameter_ForO();

    // get frame_ForO's DH matrix
    // frame_ForC's DH matrix will change during the procedure, because angle before C can rotate
    for (int i = 0; i < Population; i++){
        int numAA = particle[i].numAA;
        for (int j = 0; j < numAA; j++){
            updateDHmatrix(frame_ForO[i][j]);
        }
    }
}

void createFrameXYZ_ForO(){
    for (int i = 0; i < Population; i++){
        int numAA = particle[i].numAA;
        for (int j = 0; j < numAA; j++){
            frame_ForC[i][j].Z = std_VECTOR3d(frame_ForO[i][j].origin - frame_ForC[i][j].origin);
            frame_ForC[i][j].X = std_VECTOR3d(cross_VECTIR3d(frame_ForC[i][j].Z, frame[i][j * 3 + 2].Z));
            frame_ForC[i][j].Y = std_VECTOR3d(cross_VECTIR3d(frame_ForC[i][j].Z, frame_ForC[i][j].X));

            frame_ForO[i][j].Z = frame_ForC[i][j].Z;
            frame_ForO[i][j].X = frame_ForC[i][j].X;
            frame_ForO[i][j].Y = frame_ForC[i][j].Y;
        }
    }
}

/*
 * we can imagine a chain, from CA -> C -> O, which means if we want to know O,
 * we need to find out the transfer matrix of C (the atom before O) and pass it to O
 */
void setDH_parameter_ForO(){
    for (int i = 0; i < Population; i++){
        int numAA = particle[i].numAA;
        for (int j = 0; j < numAA; j++){
            COORDINATES *FO  = &frame_ForO[i][j];
            COORDINATES *FC  = &frame_ForC[i][j];
            COORDINATES *FCA = &frame[i][j * 3 + 2];                // bug fixed: j written wrongly to numAA

            FC->d           = dist_POINT3d(FC->origin, FCA->origin);
            FC->a           = 0;
            FC->init_theta  = findAngle(FCA->X, FC->X, FCA->Z);
            FC->alpha       = findAngle(FCA->Z, FC->Z, FC->X);

            FO->d           = dist_POINT3d(FC->origin, FO->origin);
            FO->a           = 0;
            FO->init_theta  = findAngle(FC->X, FO->X, FC->Z);
            FO->alpha       = findAngle(FC->Z, FO->Z, FO->X);
/*
            printf("FCA.X : (%lf, %lf, %lf)\n", FCA->X.x, FCA->X.y, FCA->X.z);
            printf("FCA.Z : (%lf, %lf, %lf)\n", FCA->Z.x, FCA->Z.y, FCA->Z.z);
            printf("frame for C, init theta : %lf  alpha : %lf\n\n", frame_ForC[i][j].init_theta, frame_ForC[i][j].alpha);
            */
        }
    }
}

void updateDHmatrix(COORDINATES &frameI) {
    double theta = frameI.theta + frameI.init_theta;
    double sTheta = sin(theta), cTheta = cos(theta);
    double sAlpha = sin(frameI.alpha), cAlpha = cos(frameI.alpha);
    //printf("sTheta : %lf  cTheta : %lf\n", sTheta, cTheta);
    //printf("sAlpha : %lf  cAlpha : %lf\n\n", sAlpha, cAlpha);

    double M[4][4];
    M[0][0] = cTheta;   M[0][1] = -sTheta * cAlpha;     M[0][2] = sTheta * sAlpha;      M[0][3] = frameI.a * cTheta;
    M[1][0] = sTheta;   M[1][1] = cTheta * cAlpha;      M[1][2] = -cTheta * sAlpha;     M[1][3] = frameI.a * sTheta;
    M[2][0] = 0;        M[2][1] = sAlpha;               M[2][2] = cAlpha;               M[2][3] = frameI.d;
    M[3][0] = 0;        M[3][1] = 0;                    M[3][2] = 0;                    M[3][3] = 1;


    copyMatrix(M, frameI.DH_matrix);
}

void updateParticleCoordinate(Particle &particle, int atomID, double transferMat[4][4]){
    VECTOR3d B(0, 0, 0);
    VECTOR3d ans = transferMat * B;

    particle.origin[atomID * 3 + 0] = ans.x;
    particle.origin[atomID * 3 + 1] = ans.y;
    particle.origin[atomID * 3 + 2] = ans.z;

}

/*
 *  frame[0] refers to the world coordinate in DH method, so frame[1] refers to atomID = 0
 *  for the last atom C, the key before it is not in the phi.txt too
 *  compared to original [get R transform] method, the direction of rotation angle is opposite in DH method
 *  so we need to return -theta to get the position's change
 */
double getTheta(const Particle & particle, const int &atomID){
    double theta;
    if (atomID <= 1 || atomID == particle.numAtom - 1)
        theta = 0;
    else        // atomID = 2, the first C, before it there's the first able-rotate key
        theta = (particle.Position[atomID - 2] - particle.old_position[atomID - 2]) / 180.0 * M_PI;
    return -theta;
}

VECTOR3d operator-(const POINT3d &A, const POINT3d &B){
    VECTOR3d vec;
    vec.x = A.x - B.x;
    vec.y = A.y - B.y;
    vec.z = A.z - B.z;

    return vec;
}

VECTOR3d operator-(const VECTOR3d &A, const VECTOR3d &B){
    VECTOR3d ans(A.x - B.x, A.y - B.y, A.z - B.z);
    return ans;
}

VECTOR3d operator*(double A[4][4], const VECTOR3d &B){
    double tmpVec[4] = {B.x, B.y, B.z, 1};
    double ansVec[4] = {0, 0, 0, 0};
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            ansVec[i] += A[i][j] * tmpVec[j];

    VECTOR3d ans(ansVec[0], ansVec[1], ansVec[2]);
    return ans ;
}


VECTOR3d cross_VECTIR3d(const VECTOR3d &A, const VECTOR3d &B){
    VECTOR3d ans;
    ans.x = A.y * B.z - A.z * B.y;
    ans.y = A.z * B.x - A.x * B.z;
    ans.z = A.x * B.y - A.y * B.x;

    if (VecLength(ans) == 0)        // if A // B, then just give a vector of (1, 0, 0)
        ans = VECTOR3d(1, 0, 0);

    return ans;
}

double dist_POINT3d(const POINT3d &A, const POINT3d &B){
    double x = A.x - B.x;
    double y = A.y - B.y;
    double z = A.z - B.z;
    return sqrt(x * x + y * y + z * z);
}

double findAngle(const VECTOR3d &startVec, const VECTOR3d &endVec, const VECTOR3d &refVec){
    double a = VecLength(startVec);
    double b = VecLength(endVec);
    double c = VecLength(endVec - startVec);
    //printf("a = %lf  b = %lf  c = %lf\n", a, b, c);
    double cosAngle = (a * a + b * b - c * c) / (2 * a * b);
    double angle = acos(cosAngle);
    //  judge if ans = angle or ans = -angle
    VECTOR3d direction = cross_VECTIR3d(startVec, endVec);
    double refValue = direction.x * refVec.x + direction.y * refVec.y + direction.z * refVec.z;
    isEqualToZero(refValue);
    if (refValue >= 0)
        return angle;
    else
        return -angle;
}

double VecLength(const VECTOR3d &V){
    return sqrt(V.x * V.x + V.y * V.y + V.z * V.z);
}

void copyMatrix(double src[4][4], double dst[4][4]){
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            dst[i][j] = src[i][j];
}

void isEqualToZero(double &x){
    if (fabs(x) < 0.00001)
        x = 0;
}

VECTOR3d std_VECTOR3d(const VECTOR3d &V){
    double len = VecLength(V);
    return VECTOR3d(V.x / len, V.y / len, V.z / len);
}
