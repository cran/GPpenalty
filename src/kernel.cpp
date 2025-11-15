#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;
using namespace std;


// [[Rcpp::export]]
NumericMatrix kernel_exp(const NumericMatrix& x1, const NumericMatrix& x2, const NumericVector& theta, double g) {
  int n1 = x1.nrow();
  int n2 = x2.nrow();
  int d = x1.ncol();
  // matrix k
  NumericMatrix k(n1, n2);

  for(int i=0; i<n1; i++) {
    for(int j=0; j<n2; j++) {
      double dist = 0;
      for(int l=0; l<d; l++) {
        double diff = x1(i, l) - x2(j, l);
        dist += theta[l] * diff * diff;
      }
      k(i, j) = exp(-dist);
    }
  }

  // add g to diagonal
  if(n1==n2) {
    for(int i=0; i<n1; i++) {
      k(i, i) = k(i, i) + g;
    }
  }
  return k;
}

// [[Rcpp::export]]
NumericMatrix eucli_dist(const NumericMatrix& x1, const NumericMatrix& x2) {
  int n1 = x1.nrow();
  int n2 = x2.nrow();
  int d = x1.ncol();
  NumericMatrix dist_mat(n1, n2);

  for(int i=0; i<n1; i++) {
    for(int j=0; j<n2; j++) {
      double sum_sq=0;
      for(int l=0; l<d; l++) {
        double diff = x1(i, l) - x2(j, l);
        sum_sq += diff * diff;
      }
      dist_mat(i, j) = sum_sq;
    }
  }
  return dist_mat;
}
