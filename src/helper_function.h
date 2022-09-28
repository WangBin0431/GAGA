#ifndef HELPER_FUNCTION_H
#define HELPER_FUNCTION_H

#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <stdarg.h>
#include <stdio.h>
using namespace std;

#define dsign(a,b) (b<0? -abs(a) : abs(a))

double myfmin(double ax, double bx, double(*f)(double, std::vector<Eigen::MatrixXd const*> const &), int maxitrNum, double tol, std::vector<Eigen::MatrixXd const*> const &Plist);

#endif
