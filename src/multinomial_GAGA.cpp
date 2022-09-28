//[[Rcpp::depends(RcppEigen)]]



#include <stdio.h>
#include <Rcpp.h>
#include <RcppEigen.h>
#include "helper_function.h"

using namespace std;

double Func_u_multinomial(Eigen::MatrixXd const &u, Eigen::MatrixXd const &X, Eigen::MatrixXd const &y, Eigen::MatrixXd const &b) {
  double eps = 2.2204e-16;
  int C = u.cols();
  Eigen::MatrixXd Xu = X * u;
  double v1 = -(y.leftCols(C).array() * Xu.array()).sum();
  double v2 = (Xu.array().exp().rowwise().sum() + 1).log().sum();
  double v3 = 0.5* (b.array()*(u.array().pow(2))).sum();
  return  v1 + v2 + v3;
}

double Func_lambda_multinomial(double lambda, std::vector<Eigen::MatrixXd const*> const &Plist) {
  if (Plist.size() != 5) {
    //std::cerr << "Func_lambda_logistic need 5 input parameters!" << std::endl;
    return 10000;
  }
  Eigen::MatrixXd const*u = Plist[0];
  Eigen::MatrixXd const*X = Plist[1];
  Eigen::MatrixXd const*y = Plist[2];
  Eigen::MatrixXd const*b = Plist[3];
  Eigen::MatrixXd const*d = Plist[4];

  return Func_u_multinomial((*u).array() + lambda*(*d).array(), *X, *y, *b);
}

Eigen::MatrixXd Dfu_multinomial(Eigen::MatrixXd const &u, Eigen::MatrixXd const &X,
                                Eigen::MatrixXd const &y, Eigen::MatrixXd const &b) {

  Eigen::MatrixXd Xu = X * u;
  int C = u.cols();
  Eigen::MatrixXd tmp1, tmp2;
  tmp1 = Xu.array().exp();
  tmp2 = 1 / (tmp1.array().rowwise().sum() + 1);
  return X.transpose() * (-y.leftCols(C) + (tmp1.array().colwise()*tmp2.col(0).array()).matrix()) + (u.array()*b.array()).matrix();
}

Eigen::MatrixXd getDDfu_multinomial(Eigen::MatrixXd const &u, Eigen::MatrixXd const &X,
                                    Eigen::MatrixXd const &y, Eigen::MatrixXd const &b, bool fdiag) {

  Eigen::MatrixXd Xu = X * u;
  int C = u.cols();
  int P = u.rows();
  Eigen::MatrixXd Xt = X.transpose();
  Eigen::MatrixXd invgg = Eigen::MatrixXd::Zero(P*C,P*C);
  Eigen::MatrixXd gg = Eigen::MatrixXd::Zero(P*C,1);
  Eigen::MatrixXd tmp2 = 1 + Xu.array().exp().rowwise().sum();
  Eigen::MatrixXd tmp22 = tmp2.array().pow(2);
  Eigen::MatrixXd tmp1, tmp3, tmp4, tmp5, tmp6;
  for (int ii = 0; ii < C; ii++) {
    int indexii = ii*P;
    for (int jj = 0; jj < C; jj++) {
      int indexjj = jj*P;
      if (ii == jj) {
        tmp1 = Xu.col(jj).array().exp();
        tmp4 = tmp1.array() * (tmp2.array() - tmp1.array()) / tmp22.array();
        if (fdiag == false) {
          tmp5 = Xt * (X.array().colwise() * tmp4.array().col(0)).matrix();
          tmp5.diagonal() += b.col(jj);
          invgg.block(indexjj, indexii, P, P) = tmp5;
        }
        else {
          tmp6 = (X.array().pow(2).colwise() * tmp4.array().col(0)).colwise().sum();
          gg.block(indexjj, 0, P, 1) = tmp6.transpose() + b.col(jj);
        }
      }
      else {
        if (fdiag == false) {
          tmp1 = (Xu.col(jj) + Xu.col(ii)).array().exp();
          tmp3 = Xt * (X.array().colwise()*(tmp1.array() / tmp22.array()).col(0)).matrix();
          invgg.block(indexjj, indexii, P, P) = tmp3;
          invgg.block(indexii, indexjj, P, P) = tmp3;
        }
      }
    }//for(jj in 0:C-1)
  }//for(ii in 0:C-1)
  if (fdiag == true) {
    gg = 1 / (gg.array());
    return Eigen::MatrixXd(gg.asDiagonal());
  }
  else {
    return invgg.inverse();
  }
}

Eigen::MatrixXd getEb_multinomial(Eigen::MatrixXd const &X, Eigen::MatrixXd const &y, Eigen::MatrixXd const &b, Eigen::MatrixXd const &beta, Eigen::MatrixXd &D0, int maxItr, bool fdiag) {
  int N = X.rows();
  int P = X.cols();
  int C = beta.cols();
  Eigen::MatrixXd D = D0;
  Eigen::MatrixXd u = beta;
  Eigen::MatrixXd g = Dfu_multinomial(u, X, y, b);
  Eigen::MatrixXd d;
  std::vector<Eigen::MatrixXd const*> Plist(5);
  Plist[0] = &u;
  Plist[1] = &X;
  Plist[2] = &y;
  Plist[3] = &b;
  Plist[4] = &d;

  for (int ii = 0; ii < maxItr; ii++) {
    g.resize(P*C, 1);
    d = -D * g;
    d.resize(P, C);
    double LL = myfmin(0, 2, Func_lambda_multinomial, 20, 1e-19, Plist);

    u += LL*d;
    D = getDDfu_multinomial(u, X, y, b, fdiag);

    g = Dfu_multinomial(u, X, y, b);
    if (sqrt(g.array().pow(2).sum()) / P < 1e-16) {
      break;
    }
  }
  D0 = D;
  return u;
}

// [[Rcpp::export]]
Rcpp::List cpp_multinomial_gaga(Eigen::MatrixXd X, Eigen::MatrixXd y, SEXP s_alpha, SEXP s_itrNum, SEXP s_thresh,
                                     SEXP s_flag, SEXP s_lamda_0, SEXP s_fdiag) {

  double alpha = Rcpp::as<double>(s_alpha);
  int itrNum = Rcpp::as<int>(s_itrNum);
  double thresh = Rcpp::as<double>(s_thresh);
  bool flag = Rcpp::as<bool>(s_flag);
  double lamda_0 = Rcpp::as<double>(s_lamda_0);
  bool fdiag = Rcpp::as<bool>(s_fdiag);

  bool exitflag = false;
  double eps = 1.e-19;
  int N = X.rows();
  int P = X.cols();

  int K = y.cols();
  int C = K - 1;

  Eigen::MatrixXd b, b_old, db, beta, beta_old, cov_beta, E_pow_beta, cov0;
  b = Eigen::MatrixXd::Ones(P, C).array()*lamda_0;
  b_old = b;

  int index = 1;
  for (index = 1; index <= itrNum; index++) {
    if (index == itrNum || exitflag) {
      db = b - b_old;
      b = b / alpha;
    }

    if (index == 1) {
      beta = Eigen::MatrixXd::Zero(P, C);
      cov_beta = getDDfu_multinomial(beta, X, y, b, fdiag);
    }
    int maxItr = 20;
    beta = getEb_multinomial(X, y, b, beta, cov_beta, maxItr, fdiag);
    Eigen::MatrixXd beta2 = beta.array().pow(2);
    beta2.resize(P*C,1);
    E_pow_beta = cov_beta.diagonal().array() + beta2.array();
    b = alpha / E_pow_beta.array();
    b.resize(P, C);

    if (flag && (index == itrNum || exitflag)) {
      cov0 = getDDfu_multinomial(beta, X, y, Eigen::MatrixXd::Zero(P,C), fdiag);
      Eigen::MatrixXd diagcov0 = cov0.diagonal();
      for (int k = 0; k < diagcov0.size(); k++) {
        if (E_pow_beta(k) < diagcov0(k)) beta(k) = 0;
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
