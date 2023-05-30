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
  double norm;
  double maximum;
  int kerni;
  Rcpp::NumericVector alpha;
  int type;
  double (ShotnoiseField::*kernel)(double );
  double (ShotnoiseField::*valueP)(double, double);
  double (ShotnoiseField::*valueP3)(double, double, double);
  
public:
  ShotnoiseField(SEXP snfield);
  virtual ~ShotnoiseField();
  
  double getValue(double, double);
  double getValue(double, double, double);
  double getValueSum(double, double);
  double getValueSum3(double, double, double);
  double getValueProd(double, double);
  double getValueProd3(double, double, double);
  double dist(int *, double, double);
  double dist(int *, double, double, double);
  double kernel_gauss(double r);
  double kernel_step(double r);
  double getMax();
  double getMark(int*);
};

#endif  