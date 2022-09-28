//[[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector>

using namespace Rcpp;
using namespace std;
using namespace Eigen;

Eigen::MatrixXd cpp_omp(Eigen::MatrixXd const &X, Eigen::MatrixXd const &y, int L, double eps) {
  
  int n = X.rows();
  int K = X.cols();
  Eigen::MatrixXd Xt = X.transpose();
  Eigen::MatrixXd residual = y;
  Eigen::MatrixXd proj;
  Eigen::MatrixXd X_sub;
  Eigen::MatrixXd a;
  
  std::vector<int> indx(L, 0);
  
  double tmp1;
  int pos, tmp2;
  int j = 0;
  for (j = 0; j<L; j++) {
    proj = Xt * residual;
    tmp1 = proj.array().abs().maxCoeff(&pos, &tmp2);
    indx[j] = pos;
    Xt.row(pos).setZero();
    X_sub.resize(n, j + 1);
    for (int m = 0; m<j + 1; m++) {
      X_sub.col(m) = X.col(indx[m]);
    }//for(int m=0;m<j+1;m++)
    a = (X_sub.transpose() * X_sub).ldlt().solve(X_sub.transpose() * y);
    residual = y - X_sub * a;		
    if (residual.array().pow(2).sum() <= eps) {
      break;
    }
  }//for(int j=0;j<L;j++)
  
  Eigen::MatrixXd temp = MatrixXd::Zero(K, 1);
  for (int m = 0; m<std::min(j + 1, L); m++) {
    temp(indx[m],0) = a(m, 0);
  }
  
  return temp;
}

Eigen::MatrixXd getEb_LM(Eigen::MatrixXd const &XtX, Eigen::MatrixXd const &Xty ,double sigma2, Eigen::MatrixXd const &B) {
  Eigen::MatrixXd tmp = XtX.array() + B.array()*sigma2;
  return tmp.ldlt().solve(Xty); 
}

Eigen::MatrixXd getEb_LM(Eigen::MatrixXd const &XtX, Eigen::MatrixXd const &Xty, double lambda) {
  MatrixXd tmp = XtX.array();
  tmp.diagonal().array() += lambda;
  return tmp.ldlt().solve(Xty);
}

Eigen::MatrixXd getEbb_LM(Eigen::MatrixXd const &XtX, double sigma2, Eigen::MatrixXd const &B) {
  Eigen::MatrixXd tmp = XtX.array() + B.array()*sigma2;
  tmp = tmp.inverse().array()*sigma2;
  return tmp;
}

Eigen::MatrixXd getEbb_LM(Eigen::MatrixXd const &XtX, double sigma2) {
  Eigen::MatrixXd tmp = XtX.array();
  tmp = tmp.inverse().array()*sigma2;
  return tmp;
}

Eigen::MatrixXd getEb_orth_LM(Eigen::MatrixXd const &kk, Eigen::MatrixXd const &Kty, double sigma2, Eigen::MatrixXd const &b) {
  return Kty.array() / (kk.array()+sigma2*b.array());	
}

Eigen::MatrixXd getEbb_orth_LM(Eigen::MatrixXd const &kk, double sigma2, Eigen::MatrixXd const &b) {
  return 1/(kk.array()/sigma2 + b.array());
}

Eigen::MatrixXd getEbb_orth_LM(Eigen::MatrixXd const &kk, double sigma2) {
  return 1 / (kk.array() / sigma2);
}
// [[Rcpp::export]]
Rcpp::List rcpp_lm_gaga(Eigen::MatrixXd X, Eigen::MatrixXd y,SEXP s_alpha, SEXP s_itrNum, SEXP s_thresh, SEXP s_QR_flag, 
                  SEXP s_flag, SEXP s_lamda_0, SEXP s_fix_sigma, SEXP s_sigm2_0, SEXP s_fdiag){
 
  double alpha = Rcpp::as<double>(s_alpha);
  int itrNum = Rcpp::as<int>(s_itrNum);
  double thresh = Rcpp::as<double>(s_thresh);
  bool QR_flag = Rcpp::as<bool>(s_QR_flag);
  bool flag = Rcpp::as<bool>(s_flag);
  double lamda_0 = Rcpp::as<double>(s_lamda_0);
  bool fix_sigma = Rcpp::as<bool>(s_fix_sigma);
  double sigm2_0 = Rcpp::as<double>(s_sigm2_0);
  bool fdiag = Rcpp::as<bool>(s_fdiag);
  
  bool exitflag = false;	
  double eps = 1.e-19;
  double sigm2 = 1.0;
  Eigen::MatrixXd B;
  Eigen::MatrixXd B_old;
  Eigen::VectorXd dB;
  
  if (!QR_flag) {
    int N = X.rows();
    int P = X.cols();
    sigm2 = fix_sigma ? sigm2_0 : 1.0;
    B = Eigen::MatrixXd::Identity(P, P)*lamda_0;
    B_old = B;
    
    Eigen::MatrixXd Xty = X.transpose() * y;
    Eigen::MatrixXd XtX = X.transpose() * X;
    Eigen::MatrixXd yty = y.transpose() * y;
    Eigen::MatrixXd beta, beta_old, cov_beta, E_pow_beta_M, E_pow_beta, b;
    
    for (int index = 1; index <= itrNum; index++) {
      if (sigm2 == 0)
        sigm2 = eps;
      
      if ((index == itrNum ) || exitflag) {
        dB = (B - B_old).diagonal();
        B /= alpha;
      }//if ((index == itrNum - 1) || exitflag)
      beta = getEb_LM(XtX, Xty, sigm2, B);
      cov_beta = getEbb_LM(XtX, sigm2, B);
      E_pow_beta_M = cov_beta.array() + (beta * beta.transpose()).array();
      E_pow_beta = E_pow_beta_M.diagonal();
      b = alpha / E_pow_beta.array();
      
      if (flag && (index == itrNum || exitflag)) {
        int tmpQ = 0;
        for (int k = 0; k < dB.size(); k++) {
          if (dB(k) <= 100) tmpQ += 1;
          else beta(k) = 0;
        }
        if (tmpQ == 0) {
          beta.setZero();
          return Rcpp::List::create(Rcpp::Named("itrNum") = index,
                                    Rcpp::Named("beta") = beta);
        }
        //if (tmpQ <= N) {					
        // X_sub = X[,dB<=100]; beta_sub = beta[dB<=100]; E_pow_beta_sub = E_pow_beta[dB<=100];
        Eigen::MatrixXd X_sub(N, tmpQ);
        Eigen::MatrixXd beta_sub(tmpQ, 1);
        Eigen::MatrixXd E_pow_beta_sub(tmpQ, 1);
        Eigen::VectorXi tmpidx(tmpQ);
        int tmpi = 0;
        for (int k = 0; k < dB.size(); k++) {
          if (dB(k) <= 100) {
            X_sub.col(tmpi) = X.col(k);
            beta_sub(tmpi) = beta(k);
            E_pow_beta_sub(tmpi) = E_pow_beta(k);
            tmpidx(tmpi) = k;
            tmpi += 1;
          }
        }//X_sub = X[,dB<=100]; beta_sub = beta[dB<=100]; E_pow_beta_sub = E_pow_beta[dB<=100];
        Eigen::MatrixXd CRB_sub = getEbb_LM(X_sub.transpose()*X_sub, sigm2);
        for (int k = 0; k < tmpQ; k++) {
          if (E_pow_beta_sub(k) < CRB_sub(k))beta_sub(k) = 0;
          beta(tmpidx(k)) = beta_sub(k);
        }
        //}//if (tmpQ <= N)				
        return Rcpp::List::create(Rcpp::Named("itrNum") = index,
                                  Rcpp::Named("beta") = beta);
      }//if (flag && (index == itrNum || exitflag))
      
      B_old = B;
      B.diagonal() = b;
      
      if (fix_sigma) {
        sigm2 = sigm2_0;
      }
      else {
        //sigm2 = (yty-2*t(beta)%*%Xty)/N + sum(E_pow_beta_M*XtX)/N
        sigm2 = (yty - 2 * beta.transpose() * Xty).sum() / N + (E_pow_beta_M.array()*XtX.array()).sum() / N;				
      }
      
      if (index == 1) {
        beta_old = beta;
      }
      else {
        if ((beta - beta_old).array().abs().maxCoeff() < thresh)exitflag = true;
        beta_old = beta;
      }
      
    }//for (int index = 0; index < itrNum; index++)
  }//if(!QR_flag)
  else {
    Eigen::MatrixXd beta, b, b_old, kk, Kty, EW, covW, EWW;
    Eigen::MatrixXd yty = y.transpose() * y;
    int N = X.rows();
    int P = X.cols();		
    std::vector<size_t> idx1(P);
    std::iota(idx1.begin(), idx1.end(), 0);
    if (P > N) {
      beta = cpp_omp(X, y, N, 0);
      std::sort(idx1.begin(), idx1.end(),
                [&beta](size_t index_1, size_t index_2) { return beta(index_1)*beta(index_1) > beta(index_2)*beta(index_2); });
      Eigen::MatrixXd X_sub(N,N);
      for (int k = 0; k < N; k++) {
        X_sub.col(k) = X.col(idx1[k]);
      }
      X = X_sub;
    }
    
    beta = getEb_LM(X.transpose()*X, X.transpose()*y, 100*lamda_0);
    std::vector<size_t> idx(beta.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(),
              [&beta](size_t index_1, size_t index_2) { return beta(index_1)*beta(index_1) > beta(index_2)*beta(index_2); });

    Eigen::MatrixXd sortX(X.rows(), X.cols());
    for (int k = 0; k < X.cols(); k++) {
      sortX.col(k) = X.col(idx[k]);
    }
    
    //QR 
    Eigen::HouseholderQR<Eigen::MatrixXd> qr;
    qr.compute(sortX);
    Eigen::MatrixXd K;
    Eigen::MatrixXd R;
    if (sortX.rows() > sortX.cols()) {
      Eigen::MatrixXd tmpK = qr.householderQ();
      Eigen::MatrixXd tmpR = qr.matrixQR().triangularView<Eigen::Upper>();
      K = tmpK.leftCols(sortX.cols());
      R = tmpR.topRows(sortX.cols());
    }
    else {
      K = qr.householderQ();
      R = qr.matrixQR().triangularView<Eigen::Upper>();
    }

    Eigen::MatrixXd tmp(R.diagonal().asDiagonal());
    Eigen::MatrixXd invtmp = tmp;
    invtmp.diagonal() = 1 / (invtmp.diagonal().array());
    R = invtmp * R;
    K = K * tmp;
    /*cout << "sortX:\n" << sortX << endl;
     cout << "K:\n" << K << endl;
     cout << "R:\n" << R << endl;
     cout << "||K*R - sortX||_2:\n" << (K*R - sortX).array().pow(2).sum() << endl;*/

    int N_K = K.rows();
    int P_K = K.cols();
    int N_R = R.rows();
    int P_R = R.cols();

    sigm2 = fix_sigma ? sigm2_0 : 1;
    b = lamda_0 * Eigen::MatrixXd::Ones(P_K, 1);
    b_old = b;

    kk = (K.transpose()*K).diagonal();
    Kty = K.transpose()*y;

    for (int index = 1; index <= itrNum; index++) {
      if (sigm2 == 0)
        sigm2 = eps;
      if (index == itrNum) {
        dB = b - b_old;
        b = b.array() / alpha;
      }
      EW = getEb_orth_LM(kk, Kty, sigm2, b);
      covW = getEbb_orth_LM(kk, sigm2, b);
      EWW = covW.array() + EW.array().pow(2);
      //cout << "EW:\n" << EW.transpose() << endl;
      if (flag && (index == itrNum)) {
        int tmpQ = 0;
        for (int k = 0; k < dB.size(); k++) {
          if (dB(k) <= 100) tmpQ += 1;
          else EW(k) = 0;
        }

        Eigen::MatrixXd K_sub(N_K, tmpQ);
        Eigen::MatrixXd EW_sub(tmpQ, 1);
        Eigen::MatrixXd CRB_sub(tmpQ, 1);
        Eigen::MatrixXd EWW_sub(tmpQ, 1);
        Eigen::VectorXi tmpidx(tmpQ);

        int tmpi = 0;
        for (int k = 0; k < dB.size(); k++) {
          if (dB(k) <= 100) {
            K_sub.col(tmpi) = K.col(k);
            EW_sub(tmpi) = EW(k);
            EWW_sub(tmpi) = EWW(k);
            tmpidx(tmpi) = k;
            tmpi += 1;
          }
        }//K_sub = K[,dB<=100]; EW_sub = EW[dB<=100];
        CRB_sub = getEbb_orth_LM((K_sub.transpose()*K_sub).diagonal(), sigm2);
        for (int k = 0; k < tmpQ; k++) {
          if (EWW_sub(k) < CRB_sub(k))EW_sub(k) = 0;
          EW(tmpidx(k)) = EW_sub(k);
        }//EW_sub[EWW_sub<CRB_sub] = 0; EW[dB<=100] = EW_sub;
      }//if (flag && (index == itrNum))

      b_old = b;
      b = alpha / EWW.array();
      if (fix_sigma) {
        sigm2 = sigm2_0;
      }
      else {
        sigm2 = (yty - 2 * Kty.transpose()*EW).sum() / N + (EWW.array()*kk.array()).sum() / N;
      }
    }//for (int index = 1; index <= itrNum; index++)

    Eigen::MatrixXd tmpbeta = R.inverse()*EW;
    beta.setZero();
    //beta[idx] = tmpbeta
    for (int k = 0; k < idx.size(); k++) {
      beta(idx[k]) = tmpbeta(k);
    }

    if (P > N) {
      Eigen::MatrixXd tmp(P, 1);
      tmp.setZero();
      tmp.topRows(N) = beta;

      beta = Eigen::MatrixXd::Zero(P,1);
      for (int k = 0; k < idx1.size(); k++) {
        beta(idx1[k]) = tmp(k);
      }
    }
    return Rcpp::List::create(Rcpp::Named("itrNum") = itrNum,
                              Rcpp::Named("beta") = beta);
  }
  
}
