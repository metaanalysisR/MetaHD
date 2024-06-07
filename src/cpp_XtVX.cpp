#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
NumericVector cpp_XtVX(List Qlist, List Alist, IntegerVector ind1, IntegerVector ind2) {
    int n = ind1.size();
    int m = Qlist.size();

    NumericVector sum_diag_vector(n);

    for (int k = 0; k < n; ++k) {
        double sum_diag = 0.0;
        for (int i = 0; i < m; ++i) {
            List A_sublist = Alist[i];
            NumericMatrix A1 = as<NumericMatrix>(as<NumericMatrix>(A_sublist[ind1[k] - 1]));
            NumericMatrix A2 = as<NumericMatrix>(as<NumericMatrix>(A_sublist[ind2[k] - 1]));
            arma::mat arma_A1 = as<arma::mat>(A1);
            arma::mat arma_A2 = as<arma::mat>(A2);
            arma::mat prod = arma_A1 * arma_A2;
            sum_diag += arma::accu(arma::diagvec(prod));
        }
        sum_diag_vector[k] = sum_diag;
    }
    return sum_diag_vector;
}
