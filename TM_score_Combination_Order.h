//
// Created by advanced on 18-3-19.
//

#ifndef VERSION_2_TM_SCORE_COMBINATION_ORDER_H
#define VERSION_2_TM_SCORE_COMBINATION_ORDER_H

struct COMBINATION_ORDER{
    int firRep;
    int secRep;
    int commas_amount;      // amount of :
};

extern COMBINATION_ORDER *combination_order;

void getCombinationOrder();
bool COMBINATION_sort(const COMBINATION_ORDER &ls, const COMBINATION_ORDER &rs);


#endif //VERSION_2_TM_SCORE_COMBINATION_ORDER_H
