#include <algorithm>
#include <iterator>
#include <iostream>
#include <vector>
#include <list>
#include <numeric>
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
Rcpp::NumericVector TCR_EM_counts2(Rcpp::NumericVector unique_counts, Rcpp::NumericVector counts_old, Rcpp::List t_indices, double thresh, int max_iters) {
    bool working = true;
    int iters = 0;
    Rcpp::NumericVector counts;
    while (working) {
        iters++;
        counts = clone(unique_counts);
        int listLength=t_indices.size();
        for (int i=0;i<listLength;i++) {
            Rcpp::IntegerVector idx=t_indices[i];
            int idx_len=idx.size();
            double vals[idx_len];
            for(int j=0;j<idx_len;j++){
                vals[j]=counts_old[idx[j]-1];
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
            counts_old = clone(counts);
        }
    }
    return counts;
}


// [[Rcpp::export]]
Rcpp::List make_distrs2(Rcpp::List probs_list) {
    Rcpp::List distrs=Rcpp::List::create();
    int listLength=probs_list.size();
    for (int i=0;i<listLength;i++) {
        NumericVector ps=probs_list[i];
        Rcpp::NumericVector d = {1.0};
        int vectorSize=ps.size();
        for (int k=0;k<vectorSize;k++) {
            double p=ps[k];
            int dSize=d.size();
            double y[dSize+1];
            double n[dSize+1];
            y[0]=0;
            n[dSize]=0;
            for (int j=0;j<dSize;j++) {
                y[j+1]=p*d[j];
                n[j]=(1-p)*d[j];
            }
            NumericVector tempd (dSize+1);
            for (int j = 0; j < dSize+1; j++) {
                tempd[j]=y[j]+n[j];
            }
            d=clone(tempd);
        }
        distrs.push_back(d);
    }
    return distrs;
}
