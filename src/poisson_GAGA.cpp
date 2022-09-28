//[[Rcpp::depends(RcppEigen)]]



#include <stdio.h>
#include <Rcpp.h>
#include <RcppEigen.h>
#include "helper_function.h"

using namespace std;

double Func_u_poisson(Eigen::MatrixXd const &u, Eigen::MatrixXd const &X, Eigen::MatrixXd const &y, Eigen::MatrixXd const &b) {

  Eigen::MatrixXd Xu = X * u;
  double v1 = -(y.transpose()*Xu).sum();
  double v2 = Xu.array().exp().sum();
  double v3 = 0.5* (b.array()*(u.array().pow(2))).sum();
  return  v1 + v2 + v3;
}

double Func_lambda_poisson(double lambda, std::vector<Eigen::MatrixXd const*> const &Plist) {
  if (Plist.size() != 5) {
    cerr << "Func_lambda_logistic need 5 input parameters!" << endl;
    exit(-1);
  }
  Eigen::MatrixXd const*u = Plist[0];
  Eigen::MatrixXd const*X = Plist[1];
  Eigen::MatrixXd const*y = Plist[2];
  Eigen::MatrixXd const*b = Plist[3];
  Eigen::MatrixXd const*d = Plist[4];

  return Func_u_poisson((*u).array() + lambda*(*d).array(), *X, *y, *b);
}

Eigen::MatrixXd Dfu_poisson(Eigen::MatrixXd const &u, Eigen::MatrixXd const &X,
                            Eigen::MatrixXd const &y, Eigen::MatrixXd const &b) {

  Eigen::MatrixXd Xu = X * u;
  Eigen::MatrixXd tmp = Xu.array().exp() - y.array();
  return (X.transpose() * tmp).array() + u.array()*b.array();
}

Eigen::MatrixXd getDDfu_poisson(Eigen::MatrixXd const &u, Eigen::MatrixXd const &X,
                                Eigen::MatrixXd const &y, Eigen::MatrixXd const &b, bool fdiag) {

  Eigen::MatrixXd Xu = X * u;
  if (fdiag == false) {
    Eigen::MatrixXd tmp = X.transpose()*(X.array().colwise()*Xu.array().exp().col(0)).matrix();
    tmp.diagonal() += b;
    return tmp.inverse();
  }
  else {
    Eigen::MatrixXd tmp = (X.array().pow(2).colwise()*Xu.array().exp().col(0)).colwise().sum().transpose();
    tmp += b;
    tmp = 1 / tmp.array();
    return Eigen::MatrixXd(tmp.asDiagonal());
  }
}

Eigen::MatrixXd getEb_poisson(Eigen::MatrixXd const &X, Eigen::MatrixXd const &y, Eigen::MatrixXd const &b, Eigen::MatrixXd const &beta, Eigen::MatrixXd &D0, int maxItr, bool fdiag) {
  int n = X.rows();
  int p = X.cols();
  Eigen::MatrixXd D = D0;
  Eigen::MatrixXd u = beta;
  Eigen::MatrixXd g = Dfu_poisson(u, X, y, b);
  Eigen::MatrixXd d;
  std::vector<Eigen::MatrixXd const*> Plist(5);
  Plist[0] = &u;
  Plist[1] = &X;
  Plist[2] = &y;
  Plist[3] = &b;
  Plist[4] = &d;

  for (int ii = 0; ii < maxItr; ii++) {
    d = -D * g;
    double LL = myfmin(0, 2, Func_lambda_poisson, 20, 1e-19, Plist);

    u += LL*d;
    D = getDDfu_poisson(u, X, y, b, fdiag);

    g = Dfu_poisson(u, X, y, b);
    if (sqrt(g.array().pow(2).sum()) / p < 1e-16) {
      break;
    }
  }
  D0 = D;
  return u;
}

// [[Rcpp::export]]
Rcpp::List cpp_poisson_gaga(Eigen::MatrixXd X, Eigen::MatrixXd y, SEXP s_alpha, SEXP s_itrNum, SEXP s_thresh,
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


  Eigen::MatrixXd b, b_old, db, beta, beta_old, cov_beta, E_pow_beta, cov0;
  b = Eigen::MatrixXd::Ones(p, 1).array()*lamda_0;
  b_old = b;

  int index = 1;
  for (index = 1; index <= itrNum; index++) {
    if (index == itrNum || exitflag) {
      db = b - b_old;
      b = b / alpha;
    }

    if (index == 1) {
      Eigen::MatrixXd tmpy = y;
      for (int ii = 0; ii < n; ii++) {
        if (tmpy(ii) <= 0)tmpy(ii) = 0.1;
      }
      Eigen::MatrixXd tmpM = X.transpose()*X;
      tmpM.diagonal() += b;
      beta = tmpM.ldlt().solve(X.transpose()*(tmpy.array().log().matrix()));
      //cout<<"beta:\n"<<beta.transpose()<<endl;
      cov_beta = getDDfu_poisson(beta, X, y, b, fdiag);
    }
    int maxItr = 20;
    beta = getEb_poisson(X, y, b, beta, cov_beta, maxItr, fdiag);
    E_pow_beta = cov_beta.diagonal().array() + beta.array().pow(2);
    b = alpha / E_pow_beta.array();

    if (flag && (index == itrNum || exitflag)) {
      int tmpQ = (db.array() <= 100).count();
      if (tmpQ == 0) {
        beta.setZero();
        break;
      }
      else {
        cov0 = getDDfu_poisson(beta, X, y, Eigen::MatrixXd::Zero(p, 1), fdiag);
        Eigen::MatrixXd diagcov0 = cov0.diagonal();
        for (int k = 0; k < diagcov0.size(); k++) {
          if (E_pow_beta(k) < diagcov0(k) || db(k)>20) beta(k) = 0;
        }
        break;
      }
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

