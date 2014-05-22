/*

*/

#include <vector>
#include <Rcpp.h>

#ifndef SN_H_
#define SN_H_

class ShotnoiseField {
  int n;
	Rcpp::NumericMatrix loc;
  int dim;
  double sigma;
  double max;
  int kerni;
  double (ShotnoiseField::*kernel)(double );
  
public:
  ShotnoiseField(SEXP snfield);
  virtual ~ShotnoiseField();
  
  double getValue(double, double);
  double getValue(double, double, double);
  double dist(int *, double, double);
  double dist(int *, double, double, double);
  double kernel_gauss(double r);
  double kernel_step(double r);
  double getMax();
};

#endif  