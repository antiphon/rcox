#include <Rcpp.h>
#include <Math.h>
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
  max = as<double>(field["max"]);
  
  // set the kernel
  kerni = as<int>(field["kernel"]);
  kernel = &ShotnoiseField::kernel_gauss;
  if(kerni==1) kernel = &ShotnoiseField::kernel_step;
  
  
  
  //printf("Shotnoisefield initialized [n=%i, sigma=%f].", n, sigma);
}



double ShotnoiseField::getValue(double x, double y){
  int i;
  double v=0;
  for(i=0; i < n; i++){
    v+= (this->*kernel)(dist(&i, x, y));
  }
  return v;
}


double ShotnoiseField::getValue(double x, double y, double z){
  int i;
  double v=0;
  for(i=0; i < n; i++){
    v+= (this->*kernel)(dist(&i, x, y, z));
  }
  return v;
}

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
  return max;
}
