#include "helper_function.h"
double myfmin(double ax, double bx, double(*f)(double, std::vector<Eigen::MatrixXd const*> const &), int maxitrNum, double tol, std::vector<Eigen::MatrixXd const*> const &Plist) {
  double a, b, c, d, e, eps, xm, p, q, r, tol1, tol2, u, v, w;
  double fu, fv, fw, fx, x;

  //c is the squared inverse of the golden ratio
  c = 0.5*(3. - sqrt(5.0));
  //eps is approximately the square root of the relative machine precision.
  eps = 1.0;
  goto10: eps = eps / 2.0;
  tol1 = 1.0 + eps;
  if (tol1 > 1.0) goto goto10;
  eps = sqrt(eps);


  a = ax;
  b = bx;
  v = a + c*(b - a);
  w = v;
  x = v;
  e = 0.0;
  fx = f(x, Plist);
  fv = fx;
  fw = fx;

  //main loop starts here
  int index = 0;
  goto20:
    index += 1;
  if (index == maxitrNum) goto goto90;
  xm = 0.5*(a + b);
  tol1 = eps*abs(x) + tol / 3.0;
  tol2 = 2.0*tol1;

  //check stopping criterion
  if (abs(x - xm) <= (tol2 - 0.5*(b - a))) goto goto90;

  //is golden-section necessary
  if (abs(e) <= tol1) goto goto40;

  //fit parabola
  r = (x - w)*(fx - fv);
  q = (x - v)*(fx - fw);
  p = (x - v)*q - (x - w)*r;
  q = 2.0*(q - r);
  if (q > 0.0) p = -p;
  q = abs(q);
  r = e;
  e = d;

  //is parabola acceptable
  goto30:
    if (abs(p) >= abs(0.5*q*r)) goto goto40;
    if (p <= q*(a - x)) goto goto40;
    if (p >= q*(b - x)) goto goto40;

    //a parabolic interpolation step
    d = p / q;
    u = x + d;

    //f must not be evaluated too close to ax or bx
    if ((u - a) < tol2) d = dsign(tol1, xm - x);
    if ((b - u) < tol2) d = dsign(tol1, xm - x);
    goto goto50;

    //a golden - section step
    goto40:
      if (x >= xm) e = a - x;
      if (x < xm) e = b - x;
      d = c*e;

      //f must not be evaluated too close to x
      goto50:
        if (abs(d) >= tol1) u = x + d;
        if (abs(d)< tol1) u = x + dsign(tol1, d);
        fu = f(u, Plist);

        //update  a, b, v, w, and x
        if (fu > fx) goto goto60;
        if (u >= x) a = x;
        if (u < x) b = x;
        v = w;
        fv = fw;
        w = x;
        fw = fx;
        x = u;
        fx = fu;
        goto goto20;
        goto60:
          if (u < x) a = u;
          if (u >= x) b = u;
          if (fu <= fw) goto goto70;
          if (w == x) goto goto70;
          if (fu <= fv) goto goto80;
          if (v == x) goto goto80;
          if (v == w) goto goto80;
          goto goto20;
          goto70:
            v = w;
          fv = fw;
          w = u;
          fw = fu;
          goto goto20;
          goto80:
            v = u;
          fv = fu;
          goto goto20;
          //end of main loop
          goto90:
            return x;

}
