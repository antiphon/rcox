#include <Rcpp.h>
#include <Rmath.h>
#include <stdlib.h>
#include <vector>
#include "ShotnoiseField.h"
using namespace Rcpp;
/********************************************************************************************/
  ShotnoiseField::~ShotnoiseField()
{
  }
/********************************************************************************************/
ShotnoiseField::ShotnoiseField(SEXP snfield) {
  List field(snfield);
  loc = as<NumericMatrix>(field["x"]);
  n = loc.nrow();
  dim = loc.ncol();
  sigma = as<double>(field["sigma"]);
  alpha = as<NumericVector>(field["alpha"]);
  type = as<int>(field["type"]);
  norm = 1.0; // kernel normalising constant
  // set the kernel
  kerni = as<int>(field["kernel"]);
  maximum = as<double>(field["max"]);
  // add the constants etc to have correct kernel weights
  if(kerni==1) { // box
    kernel = &ShotnoiseField::kernel_step;
    if(dim==2) norm /= (PI * sigma * sigma); // sigma = radius of disc
    else norm /= 4.0/3.0 * PI * pow(sigma, 3);
    sigma *=sigma; // to work with r^2 inputs
  }
  else{
    kernel = &ShotnoiseField::kernel_gauss;
    sigma = 2.0 * sigma * sigma;
    norm /= pow( sqrt(sigma * PI), dim);
  }
  //printf("Shotnoisefield initialized [n=%i, sigma=%f].", n, sigma);
  if(type == 0){
    valueP = &ShotnoiseField::getValueSum;
    valueP3 = &ShotnoiseField::getValueSum3;
  }
  else{
    valueP = &ShotnoiseField::getValueProd;
    valueP3 = &ShotnoiseField::getValueProd3;
  }
}

double ShotnoiseField::getValue(double x, double y){
  return (this->*valueP)(x, y);
}
double ShotnoiseField::getValue(double x, double y, double z){
  return (this->*valueP3)(x, y, z);
}
  
double ShotnoiseField::getValueSum(double x, double y){
  int i;
  double v=0;
  for(i=0; i < n; i++){
    v+= (this->*kernel)(dist(&i, x, y));
  }
  return alpha(0) * v * norm;
}
double ShotnoiseField::getValueSum3(double x, double y, double z){
  int i;
  double v=0;
  for(i=0; i < n; i++){
    v+= (this->*kernel)(dist(&i, x, y, z));
  }
  return alpha(0) * v * norm;
}

double ShotnoiseField::getValueProd(double x, double y){
  int i;
  double v=1;
  for(i=0; i < n; i++){
    v*= 1.0 + alpha(1) * (this->*kernel)(dist(&i, x, y));
  }
  return v * exp(alpha(0));
}
double ShotnoiseField::getValueProd3(double x, double y, double z){
  int i;
  double v=1;
  for(i=0; i < n; i++){
    v*= 1.0 + alpha(1) *  (this->*kernel)(dist(&i, x, y, z));
  }
  return v * exp(alpha(0));
}

  
// Note :: all distances are squared!

double ShotnoiseField::dist(int *i, double x, double y){
  double d = pow(loc(*i,0)-x, 2) + pow( loc(*i,1) - y , 2);
  return d;
}

double ShotnoiseField::dist(int *i, double x, double y, double z){
  double d = pow(loc(*i,0)-x, 2) + pow( loc(*i,1) - y , 2) + pow(loc(*i, 2) - z, 2);
  return d;
}



double ShotnoiseField::kernel_gauss(double r){ // work with squared distance
  return exp(-r/sigma);
}

double ShotnoiseField::kernel_step(double r){ // work with squared distance
  if(r < sigma) return 1;
  return 0;
}

double ShotnoiseField::getMax(){
  // the max is at mother location only for the shot-noise. Determine it:
  if(type == 0){
    double d;
    maximum = 0;
    for(int i=0; i < n; i++) {
      d = getValue(loc(i,0), loc(i,1));
      if(d > maximum) maximum = d;
    }
  }
  else{ // for the product field its at mothers if alpha(1)>0, and its exp(alpha(0)) if alpha(1)<0
    if(alpha(1) > 0) {
      double d;
      maximum = 0;
      for(int i=0; i < n; i++) {
        d = getValue(loc(i,0), loc(i,1));
        if(d > maximum) maximum = d;
      }
    }
    else{// assuming range is small enough so that mother influence doesn't cover all space
      maximum = exp(alpha(0));
    }
  }
  return maximum;
}
