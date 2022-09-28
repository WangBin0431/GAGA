//[[Rcpp::depends(RcppEigen)]]



#include <stdio.h>
#include <Rcpp.h>
#include <RcppEigen.h>
#include "helper_function.h"
//#include <numeric>
//#include <iterator>
//#include <algorithm>

using namespace std;
//using namespace Eigen;


Eigen::MatrixXd sigmod(Eigen::MatrixXd const &x) {
  return 1 / ((-x).array().exp() + 1);
}

double Func_u_logistic(Eigen::MatrixXd const &u, Eigen::MatrixXd const &X, Eigen::MatrixXd const &y, Eigen::MatrixXd const &b) {
  double eps = 2.2204e-16;
  Eigen::MatrixXd tmp = sigmod(X * u);
  Eigen::MatrixXd tmp1 = y.array() * (tmp.array() + eps).log().array();
  Eigen::MatrixXd tmp2 = (1 - y.array()).array() * (1 - tmp.array() + eps).log().array();
  return -1 * (tmp1.array() + tmp2.array()).sum() + 0.5*(b.array()*(u.array().pow(2))).sum();
}

double Func_lambda_logistic(double lambda, std::vector<Eigen::MatrixXd const*> const &Plist) {
  if (Plist.size() != 5) {
    cerr << "Func_lambda_logistic need 5 input parameters!" << endl;
    exit(-1);
  }
  Eigen::MatrixXd const*u = Plist[0];
  Eigen::MatrixXd const*X = Plist[1];
  Eigen::MatrixXd const*y = Plist[2];
  Eigen::MatrixXd const*b = Plist[3];
  Eigen::MatrixXd const*d = Plist[4];

  return Func_u_logistic((*u).array()+lambda*(*d).array(), *X, *y, *b);
}

Eigen::MatrixXd Dfu_logistic(Eigen::MatrixXd const &u, Eigen::MatrixXd const &X,
                             Eigen::MatrixXd const &y, Eigen::MatrixXd const &b) {

  Eigen::MatrixXd tmp = y.array() - sigmod(X * u).array();
  return -(X.transpose() * tmp).array() + (u.array()*b.array());
}

Eigen::MatrixXd getDDfu_logistic(Eigen::MatrixXd const &u, Eigen::MatrixXd const &X,
                                 Eigen::MatrixXd const &y, Eigen::MatrixXd const &b, bool fdiag) {


  Eigen::MatrixXd h = sigmod(X * u);
  Eigen::MatrixXd tmp = h.array()*(1 - h.array());
  Eigen::MatrixXd tmp1 = (X.array().colwise()*tmp.col(0).array()).matrix();
  //time_t start, end;
  //
  //Eigen::MatrixXd tmp1 = (X.array().colwise()*tmp.col(0).array()).matrix();
  //start = clock();
  //Eigen::MatrixXd gg = X.transpose() * tmp1;
  //end = clock();
  //std::cout << "DBSCAN time: " << (double)(end - start) / CLOCKS_PER_SEC << " s" << std::endl;         // time.h计时
  //gg.diagonal() += b;


  if (!fdiag) {
    Eigen::MatrixXd gg = X.transpose() * tmp1;
    gg.diagonal() += b;
    return gg.inverse();
  }
  else {
    //time_t start, end;
    //start = clock();
    Eigen::MatrixXd gg = (X.array()* tmp1.array()).colwise().sum();
    Eigen::MatrixXd tmp2 = 1 / (gg.transpose().array() + b.array());
    //end = clock();
    //std::cout << "DBSCAN time: " << (double)(end - start) / CLOCKS_PER_SEC << " s" << std::endl;         // time.h计时
    return Eigen::MatrixXd(tmp2.asDiagonal());
  }
}

Eigen::MatrixXd getEb_logistic(Eigen::MatrixXd const &X, Eigen::MatrixXd const &y, Eigen::MatrixXd const &b, Eigen::MatrixXd const &beta, Eigen::MatrixXd &D0, int maxItr, bool fdiag) {
  int n = X.rows();
  int p = X.cols();
  Eigen::MatrixXd D = D0;
  Eigen::MatrixXd u = beta;
  Eigen::MatrixXd g = Dfu_logistic(u, X, y, b);
  Eigen::MatrixXd d;
  std::vector<Eigen::MatrixXd const*> Plist(5);
  Plist[0] = &u;
  Plist[1] = &X;
  Plist[2] = &y;
  Plist[3] = &b;
  Plist[4] = &d;

  for (int ii = 0; ii < maxItr; ii++) {
    d = -D * g;



    double LL = myfmin(0, 2, Func_lambda_logistic, 20, 1e-19, Plist);


    u += LL*d;

    D = getDDfu_logistic(u, X, y, b, fdiag);



    g = Dfu_logistic(u, X, y, b);
    if (sqrt(g.array().pow(2).sum()) / p < 1e-16) {
      break;
    }
  }
  D0 = D;
  return u;
}

// [[Rcpp::export]]
Rcpp::List cpp_logistic_gaga(Eigen::MatrixXd X, Eigen::MatrixXd y, SEXP s_alpha, SEXP s_itrNum, SEXP s_thresh,
                                  SEXP s_flag, SEXP s_lamda_0, SEXP s_fdiag) {

  double alpha = Rcpp::as<double>(s_alpha);
  int itrNum = Rcpp::as<int>(s_itrNum);
  double thresh = Rcpp::as<double>(s_thresh);
  bool flag = Rcpp::as<bool>(s_flag);
  double lamda_0 = Rcpp::as<double>(s_lamda_0);
  bool fdiag = Rcpp::as<bool>(s_fdiag);


  bool exitflag = false;
  double eps = 1.e-19;
  int n = X.rows();
  int p = X.cols();

  Eigen::MatrixXd b, b_old, db, beta, beta_old, cov_beta, D0;

  b = Eigen::MatrixXd::Ones(p, 1).array()*lamda_0;
  b_old = b;
  int index;
  for (index = 1; index <= itrNum; index++) {
    if (index == itrNum || exitflag) {
      db = b.array() - b_old.array();
      b /= alpha;
    }

    if (index == 1) {
      beta = Eigen::MatrixXd::Zero(p, 1);
      cov_beta = getDDfu_logistic(beta, X, y, b, fdiag);
      D0 = cov_beta;
    }
    int maxItr = 20;

    //time_t start, end;
    //start = clock();

    beta = getEb_logistic(X, y, b, beta, D0, maxItr, fdiag);

    //end = clock();
    //std::cout << "DBSCAN time: " << (double)(end - start)/ CLOCKS_PER_SEC << " s" << std::endl;         // time.h计时

    cov_beta = D0;

    Eigen::MatrixXd E_pow_beta = cov_beta.diagonal().array() + beta.array().pow(2);
    b = alpha / E_pow_beta.array();

    if (flag && (index == itrNum || exitflag)) {
      Eigen::MatrixXd cov0 = getDDfu_logistic(beta, X, y, Eigen::MatrixXd::Zero(p,1), fdiag).diagonal();
      for (int k = 0; k < cov0.size(); k++) {
        if (E_pow_beta(k) < cov0(k)) beta(k) = 0;
      }
      break;
    }
    else {
      b_old = b;
    }
    if (index == 1) {
      beta_old = beta;
    }
    else {
      if ((beta - beta_old).array().abs().maxCoeff() < thresh)exitflag = true;
      beta_old = beta;
    }
  }
  return Rcpp::List::create(Rcpp::Named("itrNum") = index,
                            Rcpp::Named("beta") = beta);

}

