//
// Created by Zyj on 2018/1/21.
//

#ifndef VERSION_2_DH_ROTATION_H
#define VERSION_2_DH_ROTATION_H

#include <iostream>
#include <cmath>
#include "MOPSO.h"

class POINT3d{
public:
    POINT3d(){
        x = 0;
        y = 0;
        z = 0;
    }
    POINT3d(double x0, double y0, double z0){
        x = x0;
        y = y0;
        z = z0;
    }
public:
    double x, y, z;
};
class VECTOR3d{
public:
    VECTOR3d()
    //: x(0), y(0), z(0)
    {x = 0; y = 0; z = 0;}

    VECTOR3d(double x1, double y1, double z1)
    //: x(x1), y(y1), z(z1)
    {x = x1; y = y1; z = z1;};
public:
    double x, y, z;
};
class COORDINATES{
public:
    COORDINATES():
            theta(0), init_theta(0), d(0), a(0), alpha(0)
    {

    }
public:
    VECTOR3d X, Y, Z;
    POINT3d origin;

    // DH parameters
    double theta;
    double init_theta;
    double d;
    double a;
    double alpha;

    // DH matrix
    double DH_matrix[4][4];
};



extern COORDINATES **frame;
extern COORDINATES **frame_ForO;
extern COORDINATES **frame_ForC;

// initialize chain
void initDH_rotation();
void createChainCoordinates();
void setDHparameter();
// initialize O's frame
void initDH_paraForO();
void createFrameXYZ_ForO();
void setDH_parameter_ForO();
// update origin for particles
void convertRotationToOrigin(Particle &particle);
void getNewCoordinateForO(Particle &particle, double transferMat[4][4], const int &atomID);
// main function in DH_method
void updateDHmatrix(COORDINATES &frameI);
void updateParticleCoordinate(Particle &particle, int atomID, double transferMat[4][4]);
double getTheta(const Particle & particle, const int &atomID);

// auxiliary functions
VECTOR3d operator-(const POINT3d &A, const POINT3d &B);
VECTOR3d operator-(const VECTOR3d &A, const VECTOR3d &B);
VECTOR3d operator*(double A[4][4], const VECTOR3d &B);

VECTOR3d cross_VECTIR3d(const VECTOR3d &A, const VECTOR3d &B);
double dist_POINT3d(const POINT3d &A, const POINT3d &B);
double findAngle(const VECTOR3d &startVec, const VECTOR3d &endVec, const VECTOR3d &refVec);
double VecLength(const VECTOR3d &V);
void copyMatrix(double src[4][4], double dst[4][4]);
void isEqualToZero(double &x);
VECTOR3d std_VECTOR3d(const VECTOR3d &V);

using namespace std;


#endif //VERSION_2_DH_ROTATION_H
