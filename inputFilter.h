//
// Created by advanced on 18-3-29.
//

#ifndef VERSION_2_INPUTFILTER_H
#define VERSION_2_INPUTFILTER_H

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>//字符串分割
#include <cmath>
#include "MOPSO.h"

using namespace std;

struct RW_Data{
    int order;
};

void *getA_RwPlus_score(void *p);

/* which rep's multiply result, amount of the multiply result, the amount you need, return an array of which multiply is ok
 e.g.  repID = 1, multiAmount = 80, needAmount = 16     */
int *inputFilter(const int &repID, const int &multiAmount, const int &needAmount);

int CountLines(char *filename);
string ReadLine(char *filename,int line);
void get_RSMD_score(const int &repID, const int &producedRepAmount);
void get_RWPlus_score(const int &multiAmount);
float mean(float a[],int b);
double deviations(float x[],int y);



#endif //VERSION_2_INPUTFILTER_H
