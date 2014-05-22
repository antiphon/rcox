#include <Rcpp.h>
#include <vector>
#include "ShotnoiseField.h"
#include "helpers.h"
#include "Pp.h"

using namespace Rcpp;

// [[Rcpp::export]]
List rcox_MH(int n, NumericVector win, List snfield, int iter, int dbg) {
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
  
  // create a starting pattern
  double xnew, ynew, znew;
  for(i=0; i < n; i++){
    xnew = runif(1, win[0], win[1])(0) ;
    ynew = runif(1, win[2], win[3])(0) ;
    if(dim==3)  znew = runif(1, win[4], win[5])(0) ;
    int new_id = X.push_back(xnew, ynew, znew);
  }
  
  // then we loop
  int acc = 0;
  double E_old, E_new;
  double xold, yold, zold;
  double alpha;
  
  for(i=0; i < iter; i++) {
    j = sample_j(n); 
//    printf("potential\n");
    if(dim==2) E_old = SNF.getValue(X.getX(&j), X.getY(&j));
    else E_old = SNF.getValue(X.getX(&j), X.getY(&j), X.getZ(&j));
    xnew = runif(1, win[0], win[1])(0);
    ynew = runif(1, win[2], win[3])(0);
    xold = X.getX(&j);
    yold = X.getY(&j);
    if(dim==3) {
      znew = runif(1, win[4], win[5])(0);
      zold = X.getZ(&j);
    }
//    printf("moving\n");
    X.move_cache(&j, xnew, ynew, znew);
//    printf("potential\n");
    if(dim==2) E_new = SNF.getValue(X.getX(&j), X.getY(&j));
    else E_new = SNF.getValue(X.getX(&j), X.getY(&j), X.getZ(&j));
    
    if(E_old == 0 & E_new > 0) {alpha = 1;}
    else if(E_new == 0) { alpha = 0;}
    else {alpha = E_new/E_old; }
    
    //printf(" %f, %f, %f\n", E_old, E_new, alpha);
    
    if(runif(1)(0) < alpha) {
      acc += 1;
    }
    else {
//      printf("moving back\n");
//      X.move(&j, xold, yold, zold);
      X.move_back();
    }
    if(dbg) printf("\r %i/%i", i+1, iter);
  }
  if(dbg) printf(" MH done.\n");
  
  
  // and we are done. Compile results:
  NumericVector x(X.size()), y(X.size()), z;
  if(dim==3) z = rep(0, X.size());
  for(i=0; i < X.size(); i++) {
    x(i)=X.points.at(i).getX(); 
    y(i)=X.points.at(i).getY(); 
    if(dim==3) z(i)= X.points.at(i).getZ(); 
  }
  List xyz =   List::create(x, y, z);
  return(xyz);
}

