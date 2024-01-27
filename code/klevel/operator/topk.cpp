
#include "topk.h"
#include "oru.h"
#include "Qsim.h"
#include <vector>


void dtopk::multiple_query(level& idx, int k, int q_num, fstream& log){
    vector<vector<float>> ws; // user preference weight vectors
    generate_query(idx, q_num, ws);

    clock_t cur_time=clock();
    double visit_sum=0;
    int result_sum=0;
    for (int i=0;i<q_num;i++){
        clock_t qb=clock();
//        cout << "topk query " << i <<"("<<ws[i]<<")"<< ": " << endl;
        //log << "topk query " << i <<"("<< ws[i]<<")"<< ": " << endl;

        if (k <= 50)
            single_query(idx,k, ws[i],visit_sum, result_sum, log);

        clock_t qe=clock();
//        cout << "query time: " << (qe - qb) / (float)CLOCKS_PER_SEC << endl;
        //log << "query time: " << (qe - qb) / (float)CLOCKS_PER_SEC << endl;
    }
//    cout << "Average topk query time: " << (clock() - cur_time) / (float)CLOCKS_PER_SEC / (float) q_num << endl;
    //log << "Average kspr query time: " << (clock() - cur_time) / (float)CLOCKS_PER_SEC / (float) q_num << endl;
    log << fixed << setprecision(0) << visit_sum / q_num << "\t";
    cout << visit_sum / q_num << "\t";
    Qsim::simulate(idx, k, q_num, log, ws);
    return;
}

void dtopk::generate_query(level& idx, int q_num, vector<vector<float>>& ws_ret){
    oru::generate_query(idx, q_num, ws_ret);
}

inline bool satisfy_constr(vector<float>& point, vector<float>& constr, bool sign){
    float dot=0;
    for (int i = 0; i <point.size(); ++i) {
        dot+=point[i]*constr[i];
    }
    return !((dot>constr.back()) ^ sign);
}

bool point_in_convex(vector<float>& point, vector<halfspace> &convex){
    // check:
    // if sign:
    //     convex[:, 0:d-1] * point > convex[:, d]
    // else:
    //     convex[:, 0:d-1] * point < convex[:, d]
    if (convex.empty())
        return true;
    assert(point.size()+1==convex[0].w.size());
    for (halfspace& constr: convex) {
        if(!satisfy_constr(point, constr.w, constr.side)){
            return false;
        }
    }
    return true;
}


void dtopk::single_query(level& idx, int k, vector<float> &w, double& visit_sum, int& result_sum, fstream& log){
    kcell *cur_cell=&idx.idx[0][0];
    int i=0;
    bool flag=true;
    while(flag && i < k){
//        cout << cur_cell->curk << " " << idx.ik << endl;
        flag=false;
        if(cur_cell->curk>=idx.ik){
//            cout << "!! " << cur_cell->curk << endl;
            vector<kcell> children;
            idx.SingleCellSplit(k,*cur_cell, children);
            visit_sum += cur_cell -> r.V.size();
            for(kcell &c: children){
                idx.UpdateH(c);
                visit_sum += 1 + c.r.H.size();
                if(point_in_convex(w, c.r.H)){
                    cur_cell=&c;
                    flag=true;
                    break;
                }
            }
        }else{
            for (int & it : cur_cell->Next){
                idx.UpdateH(idx.idx[i+1][it]);
//                for (auto &w : idx.idx[i+1][it].r.H)
//                    cout << w.w[0] << " " << w.w[1] << " " << w.w[1] / w.w[0] << endl;
                visit_sum += 1 + idx.idx[i+1][it].r.H.size() + idx.idx[i+1][it].topk.size() + idx.idx[i+1][it].Stau.size();
//                visit_sum += 1 + idx.idx[i+1][it].r.H.size();
                if(point_in_convex(w, idx.idx[i+1][it].r.H)){
                    cur_cell=&idx.idx[i+1][it];
                    flag=true;
                    break;
                }
            }
        }
        ++i;
    }
    if(i<k){
//        cout<<"visit cell number: "<<visit_sum<<endl;
//        cout<<"query exists due to calculation accuracy when k="<<i<<endl;
//        cout<<"the topi options are "<<cur_cell->topk<<endl;
    }

}
