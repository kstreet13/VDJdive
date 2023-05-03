#include <algorithm>
#include <iterator>
#include <iostream>
#include <list>
#include <numeric>
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
vector<double> TCR_EM_counts2(vector<double> unique_counts, vector<double> counts_old, vector<vector<int>> t_indices, double thresh, int max_iters) {
    bool working = true;
    int iters = 0;
    vector<double> counts;
    while (working) {
        iters++;
        counts = unique_counts;
        for (vector<int> idx : t_indices) {
            int idx_len=idx.size();
            double vals[idx_len];
            for(int i=0;i<idx_len;i++){
                vals[i]=counts_old[idx[i]-1];
            }
            double s = accumulate(vals, vals+idx_len, 0.0);
            for (int ii = 0; ii < idx_len; ii++) {
                counts[idx[ii]-1] += vals[ii]/s;
            }
        }
        double diff = 0.0;
        for (int ii = 0; ii < counts.size(); ii++) {
            diff = max(diff, abs(counts[ii]-counts_old[ii]));
        }
        if ((diff < thresh) || (iters >= max_iters)) {
            working = false;
        } else {
            counts_old = counts;
        }
    }
    return counts;
}


