//
// Created by 尤肖天 on 30/12/2022.
//

#include "Qsim.h"
#include<cmath>

void Qsim::simulate(level& idx, int k, int q_num, fstream& log, vector<vector<float>>& ws){
    for (auto &w : ws){
        double tmp = 1;
        for (auto &x : w)
            tmp -= x;
        w.emplace_back(tmp);
    }
    auto item = idx.OriginD;
    int dim = item[0].size();
    long long sum = 0;

    auto calc = [&](vector<float> &p, vector<float> &w){
        double ans = 0;
        for (int i = 0; i < dim; i++)
            ans += p[i] * w[i];
        sum++;
        return ans;
    };

    //k
    //k-th element;
    sum = 0;
    vector<double> threshold;
    for (auto &w : ws){
        int l = 0;
        int r = item.size() - 1;
        while (1){
            int i = l;
            int j = r;
            double m = calc(item[i + random() % (j - i + 1)], w);
            while (i < j) {
                while (calc(item[i], w) > m) i++;
                while (calc(item[j], w) < m) j--;
                if (i <= j) {
                    swap(item[i], item[j]);
                    i++;
                    j--;
                }
            }
            if (k - 1 >= j && k - 1 < i)
                break;
            if (k - 1 < j)
                r = j;
            else
                l = i;
//            cout << l << " " << r << " " << k << endl;
        }
        double t = 1e9;
        for (int i = 0; i < k; i++)
            t = min(t, calc(item[i], w));
        threshold.emplace_back(t);
    }
    log << sum / q_num << "\t";
    cout << sum / q_num << "\t";

    //QPQ-k
    int n = item.size();
    sum = 0;
    //(double)random()/RAND_MAX
    for (int i = 0; i < q_num; i++){
        set<int, std::greater<int>> pick;
        for (int j = 0; j < k; j++){
            int x = random() % n;
            while (pick.count(x))
                x = random() % n;
            pick.emplace(x);
            sum++;
        }
        while (*pick.begin() != k - 1){
            int h = *pick.begin();
            pick.erase(h);
            bool flag = 0;
            double m = 1;
            while (!flag){
                double theta = asin(sqrt((h - k + 1) * 1.0 / n));
                int j = random() % (int)m;
                theta += j * 2 * theta;
                sum += j;
                if ((double)random()/RAND_MAX < abs(sin(theta)))
                    flag = 1;
                else
                    flag = 0;
                m *= 1.33;
                if (m > sqrt(n))
                    m = sqrt(n);
            }
            int x = random() % h;
            while (pick.count(x))
                x = random() % h;
            pick.emplace(x);
            sum++;
        }
    }

    log << sum / q_num << "\t";
    cout << sum / q_num << "\t";

    //linear
    sum = 0;
    log << idx.OriginD.size() << "\t";
    cout << idx.OriginD.size() << "\t";

    //QPQ-theta
    sum = 0;
    for (int i = 0; i < q_num; i++){
        for (int h = 1; h <= k; h++){
            bool flag = 0;
            double m = 1;
            while (!flag){
                double theta = asin(sqrt(h * 1.0 / n));
                int j = random() % ((int)m + 1);
                theta += j * 2 * theta;
                sum += j;
                double r = (double)random()/RAND_MAX;
//                cout << h << " " << j << " " << asin(h * 1.0 / n) << " " << sin(theta) << " " << r << endl;
                if (r < abs(sin(theta)))
                    flag = 1;
                else
                    flag = 0;
                m *= 1.33;
                if (m > sqrt(n))
                    m = sqrt(n);
            }
        }
    }
    log << sum / q_num << endl;
    cout << sum / q_num << endl;
}
