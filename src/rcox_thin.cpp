#include <Rcpp.h>
#include <vector>
#include "ShotnoiseField.h"
#include "helpers.h"
#include "Pp.h"

using namespace Rcpp;

// [[Rcpp::export]]
List rcox_thin(NumericVector win, List snfield, int dbg) {
  RNGScope scope;
  
  int i, j;
  
  int dim = 2;
  if(win.size()>4) dim = 3;
  
  ShotnoiseField SNF(snfield);
  
  //initial pattern
  std::vector<double> window;
  for(j=0; j < win.size(); j++) window.push_back(win(j));
  Pp X(window, 0);
  
  //if(blocking > 0) X.start_blocking(blocking);
  
  // Determine the maximum of the field for dominating Poisson
  double maximum = SNF.getMax();
  double Vol = X.getArea();
  // create a starting pattern from uniform Poisson
  int n = rpois(1, Vol * maximum)(0);
  if(dbg)Rprintf("Vol=%f, max = %f, n0=%i\n", Vol, maximum, n);
  
  double xnew, ynew, znew;
  for(i=0; i < n; i++){
    xnew = runif(1, win[0], win[1])(0) ;
    ynew = runif(1, win[2], win[3])(0) ;
    if(dim==3)  znew = runif(1, win[4], win[5])(0) ;
    X.push_back(xnew, ynew, znew);
  }
  
  // then we thin
  double p, v;
  IntegerVector keep(n);
  int total=0;
  for(i=0; i < n; i++) {
    if(dim==2) v = SNF.getValue(X.getX(&i), X.getY(&i));
    else v = SNF.getValue(X.getX(&i), X.getY(&i), X.getZ(&i));
    p = v/maximum;
    if(runif(1)[0] < p){ // keep
      keep(i) = 1;
      total++;
    }
    if(dbg) Rprintf("\r %i/%i", i+1, n);
  }
  if(dbg) Rprintf(" thinning done, kept %i.\n", total);
  
  
  // and we are done. Compile results:
  NumericVector x(total), y(total), z;
  int l=0;
  if(dim==3) z = rep(0, total);
  for(i=0; i < X.size(); i++) {
    if(keep(i)){
      x(l)=X.points.at(i).getX(); 
      y(l)=X.points.at(i).getY(); 
      if(dim==3) z(l)= X.points.at(i).getZ(); 
      l++;
    }
  }
  List xyz =   List::create(x, y, z);
  return(xyz);
}

