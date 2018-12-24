//
// Created by advanced on 17-1-16.
//

#include "MOPSO.h"

time_t starTime;
time_t endTime;

void setTime(){
    starTime = time(NULL);
    srand(time(NULL));
}

void getEndTime(){
    endTime = time(NULL);
}