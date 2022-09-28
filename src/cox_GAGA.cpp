//[[Rcpp::depends(RcppEigen)]]



#include <stdio.h>
#include <Rcpp.h>
#include <RcppEigen.h>
#include "helper_function.h"

using namespace std;

Eigen::MatrixXd mrank(Eigen::MatrixXd const& sorty) {
  Eigen::MatrixXd atrisk = Eigen::MatrixXd::Zero(sorty.size(), 1);
  atrisk(0) = 0;
  if (sorty.size() == 1)return atrisk;

  for (int ii = 1; ii < sorty.size(); ii++) {
    if (sorty(ii) == sorty(ii - 1)) {
      atrisk(ii) = atrisk(ii - 1);
    }
    else {
      atrisk(ii) = ii;
    }
  }
  return atrisk;
}

double Func_u_COX(Eigen::MatrixXd const &u, Eigen::MatrixXd const &X, Eigen::MatrixXd const &sorty, Eigen::MatrixXd const &freq, Eigen::MatrixXd const &cens, Eigen::MatrixXd const &atrisk, Eigen::MatrixXd const &b) {
  /*Eigen::MatrixXd obsfreq = freq.array() * (1 - cens.array());
   Eigen::MatrixXd Xu = X * u;
   Eigen::MatrixXd tmp = freq.array()*Xu.array().exp();
   Eigen::MatrixXd risksum = tmp;
   for (int k = 0; k < tmp.size(); k++) {
   risksum(k) = tmp(atrisk(k));
   }

   double v1 = -(obsfreq.transpose() * (Xu.array() - risksum.array().log()).matrix()).sum();
   double v2 = 0.5* (b.array()*(u.array().pow(2))).sum();
   cout << v1 + v2 << endl;
   return  v1 + v2 ;*/

  Eigen::MatrixXd obsfreq = freq.array() * (1 - cens.array());
  Eigen::MatrixXd Xu = X * u;
  Eigen::MatrixXd r = Xu.array().exp();
  Eigen::MatrixXd tmp1 = freq.array()*r.array();
  Eigen::MatrixXd tmp0 = tmp1;
  int kk = tmp0.size() - 2;
  while (kk >= 0) {
    tmp0(kk) += tmp0(kk + 1);
    kk -= 1;
  }
  Eigen::MatrixXd risksum = tmp0;
  for (int k = 0; k < tmp0.size(); k++) {
    risksum(k) = tmp0(atrisk(k));
  }

  double v1 = -(obsfreq.transpose() * (Xu.array() - risksum.array().log()).matrix()).sum();
  double v2 = 0.5* (b.array()*(u.array().pow(2))).sum();
  double v = v1 + v2;
  //cout << v << endl;
  return v;
}

double Func_lambda_COX(double lambda, std::vector<Eigen::MatrixXd const*> const &Plist) {
  if (Plist.size() != 8) {
    //std::cerr << "Func_lambda_logistic need 8 input parameters!" << std::endl;
    return 10000;
  }
  Eigen::MatrixXd const*u = Plist[0];
  Eigen::MatrixXd const*X = Plist[1];
  Eigen::MatrixXd const*sorty = Plist[2];
  Eigen::MatrixXd const*freq = Plist[3];
  Eigen::MatrixXd const*cens = Plist[4];
  Eigen::MatrixXd const*atrisk = Plist[5];
  Eigen::MatrixXd const*b = Plist[6];
  Eigen::MatrixXd const*d = Plist[7];

  return Func_u_COX((*u).array() + lambda*(*d).array(), *X, *sorty, *freq, *cens, *atrisk, *b);
}

Eigen::MatrixXd getDDfu_COX(Eigen::MatrixXd const &u, Eigen::MatrixXd const &X,
                            Eigen::MatrixXd const &sorty, Eigen::MatrixXd const &freq,
                            Eigen::MatrixXd const &cens, Eigen::MatrixXd const &atrisk,
                            Eigen::MatrixXd const &b, bool fdiag) {

  Eigen::MatrixXd gg;

  Eigen::MatrixXd obsfreq = freq.array() * (1 - cens.array());
  Eigen::MatrixXd Xu = X * u;
  Eigen::MatrixXd r = Xu.array().exp();
  Eigen::MatrixXd tmp1 = freq.array()*r.array();
  Eigen::MatrixXd tmp0 = tmp1;
  int kk = tmp0.size() - 2;
  while (kk >= 0) {
    tmp0(kk) += tmp0(kk + 1);
    kk -= 1;
  }
  Eigen::MatrixXd risksum = tmp0;
  for (int k = 0; k < tmp0.size(); k++) {
    risksum(k) = tmp0(atrisk(k));
  }

  int n = X.rows();
  int p = X.cols();

  Eigen::MatrixXd Xr = X.array().colwise() * tmp1.array().col(0);
  Eigen::MatrixXd Xrsum = Xr;

  int ii = Xr.rows() - 2;
  while (ii >= 0) {
    Xr.row(ii) += Xr.row(ii + 1);
    ii--;
  }

  for (int k = 0; k < Xr.rows(); k++) {
    Xrsum.row(k) = Xr.row(atrisk(k));
  }

  Eigen::MatrixXd A = Xrsum.array().colwise() / risksum.array().col(0);

  if (fdiag == false) {
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(n, p*(p + 1) / 2);

    /*  t1 = rep(1:p, p)
     t2 = sort(t1)
     XXr = X[, t1] * X[, t2] * (tmp1)
     XXrsum = apply(XXr, 2, revcumsum)
     XXrsum = XXrsum[atrisk, ] / risksum
     gg = t(obsfreq) % *%XXrsum
     dim(gg) = c(p, p)
     gg = ginv(gg - t(A) % *%(A*obsfreq) + diag(b))*/
    int kk = 0;
    for (int ii = 0; ii < p; ii++) {
      for (int jj = 0; jj <= ii; jj++) {
        M.col(kk++) = X.col(ii).array()*X.col(jj).array() * tmp1.array().col(0);
      }
    }
    int ii = M.rows() - 2;
    while (ii >= 0) {
      M.row(ii) += M.row(ii + 1);
      ii--;
    }
    Eigen::MatrixXd XXrsum = M;
    for (int k = 0; k < XXrsum.rows(); k++) {
      XXrsum.row(k) = M.row(atrisk(k));
    }
    XXrsum.array().colwise() /= risksum.array().col(0);
    Eigen::MatrixXd tmp = obsfreq.transpose() * XXrsum;

    gg = Eigen::MatrixXd::Zero(p, p);
    kk = 0;
    for (int ii = 0; ii < p; ii++) {
      for (int jj = 0; jj <= ii; jj++) {
        int index1 = ii*p + jj;
        int index2 = jj*p + ii;
        if (index1 != index2) {
          gg(index1) = gg(index2) = tmp(kk);
        }
        else {
          gg(index1) = tmp(kk);
        }
        kk++;
      }
    }
    gg -= A.transpose() * (A.array().colwise() * obsfreq.array().col(0)).matrix();
    gg.diagonal() += b;
    gg = gg.inverse();
  }
  else {
    /*  XXr = X*X*tmp1
     XXrsum = apply(XXr, 2, revcumsum)
     XXrsum = XXrsum[atrisk, ] / risksum
     gg = t(obsfreq) % *%XXrsum
     gg = diag(as.numeric(1 / (gg - apply(A*A*obsfreq, 2, sum) + b)))*/
    Eigen::MatrixXd XXr = X.array().pow(2).colwise()*tmp1.array().col(0);
    int ii = XXr.rows() - 2;
    while (ii >= 0) {
      XXr.row(ii) += XXr.row(ii + 1);
      ii--;
    }
    Eigen::MatrixXd XXrsum = XXr;
    for (int k = 0; k < XXrsum.rows(); k++) {
      XXrsum.row(k) = XXr.row(atrisk(k));
    }
    XXrsum.array().colwise() /= risksum.array().col(0);
    Eigen::MatrixXd tmp = obsfreq.transpose() * XXrsum;
    Eigen::MatrixXd tmp2 = (A.array().pow(2).colwise()*obsfreq.array().col(0)).colwise().sum();
    gg = Eigen::MatrixXd((1 / (tmp.transpose() - tmp2.transpose() + b).array()).matrix().asDiagonal());
  }

  return gg;
}


double negloglike_COX(Eigen::MatrixXd const &u, Eigen::MatrixXd const &X,
                      Eigen::MatrixXd const &sorty, Eigen::MatrixXd const &freq,
                      Eigen::MatrixXd const &cens, Eigen::MatrixXd const &atrisk,
                      Eigen::MatrixXd const &b, bool fdiag,
                      Eigen::MatrixXd &g, Eigen::MatrixXd &gg) {

  Eigen::MatrixXd obsfreq = freq.array() * (1 - cens.array());
  Eigen::MatrixXd Xu = X * u;
  Eigen::MatrixXd r = Xu.array().exp();
  Eigen::MatrixXd tmp1 = freq.array()*r.array();
  Eigen::MatrixXd tmp0 = tmp1;
  int kk = tmp0.size() - 2;
  while (kk >= 0) {
    tmp0(kk) += tmp0(kk + 1);
    kk -= 1;
  }
  Eigen::MatrixXd risksum = tmp0;
  for (int k = 0; k < tmp0.size(); k++) {
    risksum(k) = tmp0(atrisk(k));
  }

  double v1 = -(obsfreq.transpose() * (Xu.array() - risksum.array().log()).matrix()).sum();
  double v2 = 0.5* (b.array()*(u.array().pow(2))).sum();
  double v = v1 + v2;

  int n = X.rows();
  int p = X.cols();

  Eigen::MatrixXd Xr = X.array().colwise() * tmp1.array().col(0);
  Eigen::MatrixXd Xrsum = Xr;

  int ii = Xr.rows() - 2;
  while (ii >= 0) {
    Xr.row(ii) += Xr.row(ii + 1);
    ii--;
  }

  for (int k = 0; k < Xr.rows(); k++) {
    Xrsum.row(k) = Xr.row(atrisk(k));
  }

  Eigen::MatrixXd A = Xrsum.array().colwise() / risksum.array().col(0);

  g = (-obsfreq.transpose() * (X - A)).transpose().array() + u.array()*b.array();

  if (fdiag == false) {
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(n, p*(p + 1) / 2);

    /*  t1 = rep(1:p, p)
     t2 = sort(t1)
     XXr = X[, t1] * X[, t2] * (tmp1)
     XXrsum = apply(XXr, 2, revcumsum)
     XXrsum = XXrsum[atrisk, ] / risksum
     gg = t(obsfreq) % *%XXrsum
     dim(gg) = c(p, p)
     gg = ginv(gg - t(A) % *%(A*obsfreq) + diag(b))*/
    int kk = 0;
    for (int ii = 0; ii < p; ii++) {
      for (int jj = 0; jj <= ii; jj++) {
        M.col(kk++) = X.col(ii).array()*X.col(jj).array() * tmp1.array().col(0);
      }
    }
    int ii = M.rows() - 2;
    while (ii >= 0) {
      M.row(ii) += M.row(ii + 1);
      ii--;
    }
    Eigen::MatrixXd XXrsum = M;
    for (int k = 0; k < XXrsum.rows(); k++) {
      XXrsum.row(k) = M.row(atrisk(k));
    }
    XXrsum.array().colwise() /= risksum.array().col(0);
    Eigen::MatrixXd tmp = obsfreq.transpose() * XXrsum;

    gg = Eigen::MatrixXd::Zero(p, p);
    kk = 0;
    for (int ii = 0; ii < p; ii++) {
      for (int jj = 0; jj <= ii; jj++) {
        int index1 = ii*p + jj;
        int index2 = jj*p + ii;
        if (index1 != index2) {
          gg(index1) = gg(index2) = tmp(kk);
        }
        else {
          gg(index1) = tmp(kk);
        }
        kk++;
      }
    }
    gg -= A.transpose() * (A.array().colwise() * obsfreq.array().col(0)).matrix();
    gg.diagonal() += b;
    gg = gg.inverse();
  }
  else {
    /*  XXr = X*X*tmp1
     XXrsum = apply(XXr, 2, revcumsum)
     XXrsum = XXrsum[atrisk, ] / risksum
     gg = t(obsfreq) % *%XXrsum
     gg = diag(as.numeric(1 / (gg - apply(A*A*obsfreq, 2, sum) + b)))*/
    Eigen::MatrixXd XXr = X.array().pow(2).colwise()*tmp1.array().col(0);
    int ii = XXr.rows() - 2;
    while (ii >= 0) {
      XXr.row(ii) += XXr.row(ii + 1);
      ii--;
    }
    Eigen::MatrixXd XXrsum = XXr;
    for (int k = 0; k < XXrsum.rows(); k++) {
      XXrsum.row(k) = XXr.row(atrisk(k));
    }
    XXrsum.array().colwise() /= risksum.array().col(0);
    Eigen::MatrixXd tmp = obsfreq.transpose() * XXrsum;
    Eigen::MatrixXd tmp2 = (A.array().pow(2).colwise()*obsfreq.array().col(0)).colwise().sum();
    gg = Eigen::MatrixXd((1 / (tmp.transpose() - tmp2.transpose() + b).array()).matrix().asDiagonal());

  }

  return v;
}

Eigen::MatrixXd getEb_COX(Eigen::MatrixXd const &u0, Eigen::MatrixXd const &X,
                          Eigen::MatrixXd const &sorty, Eigen::MatrixXd const &freq,
                          Eigen::MatrixXd const &cens, Eigen::MatrixXd const &atrisk,
                          Eigen::MatrixXd const &b, Eigen::MatrixXd &D0, int maxItr, bool fdiag) {

  int n = X.rows();
  int p = X.cols();

  Eigen::MatrixXd g, D, u;
  u = u0;
  double v;
  v = negloglike_COX(u, X, sorty, freq, cens, atrisk, b, fdiag, g, D);

  Eigen::MatrixXd d;
  std::vector<Eigen::MatrixXd const*> Plist(8);
  Plist[0] = &u;
  Plist[1] = &X;
  Plist[2] = &sorty;
  Plist[3] = &freq;
  Plist[4] = &cens;
  Plist[5] = &atrisk;
  Plist[6] = &b;
  Plist[7] = &d;

  for (int ii = 0; ii < maxItr; ii++) {
    d = -D * g;

    double LL = myfmin(0, 2, Func_lambda_COX, 20, 1e-19, Plist);

    u += LL*d;
    v = negloglike_COX(u, X, sorty, freq, cens, atrisk, b, fdiag, g, D);

    if (sqrt(g.array().pow(2).sum()) / p < 1e-16) {
      break;
    }
  }
  D0 = D;
  return u;
}


// [[Rcpp::export]]
Rcpp::List cpp_COX_gaga(Eigen::MatrixXd X, Eigen::MatrixXd y, Eigen::MatrixXd cens, double alpha = 2, int itrNum = 50, double thresh = 0.001,
                             bool flag = true, double lamda_0 = 0.5, bool fdiag = true) {

  bool exitflag = false;
  double eps = 1.e-19;
  int n = X.rows();
  int p = X.cols();

  std::vector<size_t> idx(y.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(),
            [&y](size_t index_1, size_t index_2) { return y(index_1) < y(index_2); });

  Eigen::MatrixXd sortX(n, p);
  Eigen::MatrixXd sorty(n, 1);
  Eigen::MatrixXd sortcens(n, 1);
  for (int k = 0; k < X.rows(); k++) {
    sortX.row(k) = X.row(idx[k]);
    sorty(k) = y(idx[k]);
    sortcens(k) = cens(idx[k]);
  }
  //cout << sorty << endl;
  X = sortX;
  cens = sortcens;

  Eigen::MatrixXd freq = Eigen::MatrixXd::Ones(n, 1);
  double sumf = max(1.0, freq.sum());
  Eigen::MatrixXd baseX = (freq.transpose() * X).array() / sumf;
  X.array().rowwise() -= baseX.array().row(0);

  Eigen::MatrixXd atrisk = mrank(sorty);

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
      Eigen::MatrixXd stdX = ((freq.transpose() * (X.array().pow(2)).matrix()).array() / sumf).array().sqrt().transpose();
      beta = Eigen::MatrixXd::Zero(p, 1);
      for (int k = 0; k < p; k++) {
        if (stdX(k) != 0)beta(k) = 0.01 / stdX(k);
      }

    }

    int maxItr = 20;
    beta = getEb_COX(beta, X, sorty, freq, cens, atrisk, b, cov_beta, maxItr, fdiag);
    //cout<<beta<<endl;
    E_pow_beta = cov_beta.diagonal().array() + beta.array().pow(2);

    b = alpha / E_pow_beta.array();

    if (flag && (index == itrNum || exitflag)) {
      int tmpQ = (db.array() <= 100).count();
      if (tmpQ == 0) {
        beta.setZero();
        break;
      }
      else {
        cov0 = getDDfu_COX(beta, X, sorty, freq, cens, atrisk, Eigen::MatrixXd::Zero(p, 1), fdiag);

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









