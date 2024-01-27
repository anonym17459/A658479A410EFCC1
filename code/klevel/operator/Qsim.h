//
// Created by 尤肖天 on 30/12/2022.
//

#ifndef K_LEVEL_QSIM_H
#define K_LEVEL_QSIM_H
#include <vector>
#include "topk.h"
#include "oru.h"
using namespace std;


class Qsim {
public:
    static void simulate(level& idx, int k, int q_num, fstream& log, vector<vector<float>>& ws);
};


#endif //K_LEVEL_QSIM_H
