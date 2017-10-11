#include <Rcpp.h>
#include <vector>
#include "ShotnoiseField.h"

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector evaluate_lambda_c(List snfield, NumericMatrix x) {
  int i;
  int dim = x.ncol();
  int n = x.nrow();
  ShotnoiseField SNF(snfield);
  
  NumericVector v(n);
  //if(blocking > 0) X.start_blocking(blocking);
  if(dim==2)
    for(i=0; i < n; i++){
      v(i) = SNF.getValue(x(i,0), x(i,1));
    }
  if(dim==3)
    for(i=0; i < n; i++){
      v(i) = SNF.getValue(x(i,0), x(i,1), x(i,2));
    }
    
  return(v);
}

