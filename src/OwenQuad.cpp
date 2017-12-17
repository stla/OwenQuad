// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
using namespace Rcpp;
#include <boost/multiprecision/float128.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/owens_t.hpp>
#include <cmath>
#include <climits>


namespace mp = boost::multiprecision;
namespace m = boost::math;

const double one_div_root_two = m::constants::one_div_root_two<double>();
const double root_two_pi = m::constants::root_two_pi<double>();
const mp::float128 root_two_pi128 = m::constants::root_two_pi<mp::float128>();
const mp::float128 one_div_root_two_pi128 = m::constants::one_div_root_two_pi<mp::float128>();
const mp::float128 one_div_root_two128 = m::constants::one_div_root_two<mp::float128>();

double pnorm64(double q){
  if(std::isnan(q)){
    return nan("");
  }
  // if(fabs(q) > DBL_MAX){ // inutile (à confirmer mais je suis sûr)
  //   return q > 0 ? 1 : 0;
  // }
  return m::erfc(-q * one_div_root_two)/2.0;
}

mp::float128 dnorm128(mp::float128 x){
  return mp::exp(-x*x/2) * one_div_root_two_pi128;
}

mp::float128 pnorm128(mp::float128 q){
  if(fabs(q) > DBL_MAX){
    return q>0 ? mp::float128(1) : mp::float128(0);
  }
  return m::erfc(-q * one_div_root_two128)/2;
}

int sign(double x){
  return (std::signbit(x) ? -1 : 1);
}

//********* Owen T-function **************************************************//
// [[Rcpp::export]]
double owent(double h, double a){
  return m::owens_t(h,a);
}
//****************************************************************************//

// ------ Student CDF ------------------------------------------------------- //
double* studentCDF_C(double q, int nu, NumericVector delta, size_t J){
  const double a = sign(q)*sqrt(q*q/nu);
  const double sb = sqrt(nu/(nu+q*q));
  double* C = new double[J];
  size_t j;
  for(j=0; j<J; j++){
    C[j] = 2*owent(delta[j] * sb, a) + pnorm64(-delta[j]*sb);
  }
  return C;
}

// [[Rcpp::export]]
NumericVector studentCDF(double q, int nu, NumericVector delta){
  int J = delta.size();
  NumericVector out(J);
  if(nu > INT_MAX){ // nu=Inf fait planter R; dans Haskell il n'y a pas de Inf integer
    for(int j=0; j<J; j++){
      out[j] = pnorm64(q - delta[j]);
    }
    return out;
  }
  if(fabs(q)>DBL_MAX){ // faire dans R?
    for(int j=0; j<J; j++){
      out[j] = fabs(delta[j]) > DBL_MAX ?
        (std::signbit(q) == std::signbit(delta[j]) ?
          nan("") :
          (std::signbit(q) ? 0 : 1)) :
      (std::signbit(q) ? 0 : 1);
    }
    return out;
  }
  if(nu==1){
    double* C = studentCDF_C(q, nu, delta, J);
    for(int j=0; j<J; j++){
      out[j] = C[j];
    }
    delete[] C;
    return out;
  }
  const mp::float128 qq(q*q);
  const mp::float128 a = sign(q)*mp::sqrt(qq/nu);
  const mp::float128 b = nu/(nu+qq);
  const mp::float128 sb = mp::sqrt(b);
  std::vector<mp::float128> dsb(J);
  int j;
  for(j=0; j<J; j++){
    dsb[j] = delta[j] * sb;
  }
  mp::float128 M[nu-1][J];
  for(j=0; j<J ; j++){
    M[0][j] = a * sb * dnorm128(dsb[j]) * pnorm128(a*dsb[j]);
  }
  if(nu>2){
    for(j=0; j<J; j++){
      M[1][j] = b * (delta[j] * a * M[0][j] + a * dnorm128(delta[j]) * one_div_root_two_pi128);
    }
    if(nu>3){
      std::vector<mp::float128> A(nu-3);
      A[0] = 1.0;
      int k;
      if(nu>4){
        for(k=1; k<nu-3; k++){
          A[k] = 1.0/k/A[k-1];
        }
      }
      for(k=2; k<nu-1; k++){
        for(j=0; j<J; j++){
          M[k][j] = (k-1) * b * (A[k-2] * delta[j] * a * M[k-1][j] + M[k-2][j]) / k;
        }
      }
    }
  }
  if(nu%2 == 1){
    double* C = studentCDF_C(q, nu, delta, J);
    std::vector<mp::float128> sum(J);
    int i;
    for(i=1; i<nu-1; i+=2){
      for(j=0; j<J; j++){
        sum[j] += M[i][j];
      }
    }
    for(j=0; j<J; j++){
      out[j] = C[j] + 2*sum[j].convert_to<double>();
    }
    delete[] C;
    return out;
  }
  int i;
  std::vector<mp::float128> sum(J);
  for(i=0; i<nu-1; i+=2){
    for(j=0; j<J; j++){
      sum[j] += M[i][j];
    }
  }
  for(j=0; j<J; j++){
    out[j] = pnorm64(-delta[j]) + (root_two_pi128*sum[j]).convert_to<double>();
  }
  return out;
}

// ------- Owen Q-function -------------------------------------------------- //
double* OwenQ1_C(int nu, double t, NumericVector delta, NumericVector R, size_t J){
  const double a = sign(t)*sqrt(t*t/nu);
  const double b = nu/(nu+t*t);
  const double sb = sqrt(b);
  const double ab = fabs(t) > DBL_MAX ? 0 : sqrt(nu) * 1/(nu/t+t);
  double* C = new double[J];
  size_t i;
  for(i=0; i<J; i++){
    double C1 = owent(delta[i]*sb, a);
    double C2 = owent(R[i], a-delta[i]/R[i]);
    double C3 = owent(delta[i]*sb, (ab-R[i]/delta[i])/b);
    C[i] = pnorm64(R[i]) - (delta[i] >= 0) + 2*(C1 - C2 - C3);
  }
  return C;
}

// [[Rcpp::export]]
NumericVector OwenQ1(int nu, double t, NumericVector delta, NumericVector R){
  const int J = delta.size();
  NumericVector out(J);
  if(t > DBL_MAX){
    for(int j=0; j<J; j++){
      out[j] = m::gamma_p(0.5*nu, 0.5*R[j]*R[j]);
    }
    return out;
  }
  if(t < DBL_MIN || nu > INT_MAX){
    // for(int j=0; j<J; j++){
    //   out[j] = 0.0;
    // }
    return out;
  }
  if(nu == 1){
    double* C = OwenQ1_C(nu, t, delta, R, J);
    for(int j=0; j<J; j++){
      out[j] = C[j];
    }
    delete[] C;
    return out;
  }
  const mp::float128 tt(t*t);
  const mp::float128 a = sign(t)*mp::sqrt(tt/nu);
  const mp::float128 b = nu/(nu+tt);
  const mp::float128 sb = mp::sqrt(b);
  mp::float128 ab;
  mp::float128 asb;
  if(fabs(t) > DBL_MAX){ // faire un "? :"
    ab = 0;
    asb = sign(t);
  }else{
    ab = a*b;
    asb = sign(t)*mp::sqrt(tt/(nu+tt));
  }
  mp::float128 dnormdsb[J];
  mp::float128 dabminusRoversb[J];
  mp::float128 dnormR[J];
  const int n = nu-1;
  mp::float128 H[n][J];
  mp::float128 M[n][J];
  int j;
  for(j=0; j<J; j++){
    dnormdsb[j] = dnorm128(delta[j] * sb);
    dabminusRoversb[j] = (delta[j]*ab - R[j])/sb;
    dnormR[j] = dnorm128(R[j]);
    H[0][j] = -dnormR[j] * pnorm128(a*R[j]-delta[j]);
    M[0][j] = asb * dnormdsb[j] * (pnorm128(delta[j]*asb) - pnorm128(dabminusRoversb[j]));
  }
  if(nu >= 3){
    for(j=0; j<J; j++){
      H[1][j] = R[j] * H[0][j];
      M[1][j] = delta[j]*ab*M[0][j] + ab * dnormdsb[j] *
        (dnorm128(delta[j]*asb) - dnorm128(dabminusRoversb[j]));
    }
    if(nu >= 4){
      mp::float128 A[n];
      mp::float128 L[n-2][J];
      A[0] = 1;
      A[1] = 1;
      for(j=0; j<J; j++){
        L[0][j] = ab * R[j] * dnormR[j] * dnorm128(a*R[j]-delta[j])/2;
      }
      int k;
      for(k=2; k<n; k++){
        A[k] = 1.0/k/A[k-1];
      }
      if(nu >= 5){
        for(k=1; k<n-2; k++){
          for(j=0; j<J; j++){
            L[k][j] = A[k+2] * R[j] * L[k-1][j];
          }
        }
      }
      for(k=2; k<n; k++){
        for(j=0; j<J; j++){
          H[k][j] = A[k] * R[j] * H[k-1][j];
          M[k][j] = (k-1.0)/k * (A[k-2] * delta[j] * ab * M[k-1][j] + b*M[k-2][j]) - L[k-2][j];
        }
      }
    }
  }
  std::vector<mp::float128> sum(J);
  int i;
  if(nu % 2 == 0){
    for(i=0; i<n; i+=2){
      for(j=0; j<J; j++){
        sum[j] += M[i][j]+H[i][j];
      }
    }
    for(j=0; j<J; j++){
      out[j] = pnorm64(-delta[j]) + root_two_pi*sum[j].convert_to<double>();
    }
    return out;
  }else{
    for(i=1; i<n; i+=2){
      for(j=0; j<J; j++){
        sum[j] += M[i][j]+H[i][j];
      }
    }
    double* C = OwenQ1_C(nu, t, delta, R, J);
    for(j=0; j<J; j++){
      out[j] = C[j] + 2*sum[j].convert_to<double>();
    }
    delete[] C;
    return out;
  }
}

// --- Owen cumulative function 4 ------------------------------------------- //
double* OwenCDF4_C(int nu, double t1, double t2, NumericVector delta1, NumericVector delta2, size_t J){
  const double a1 = sign(t1)*sqrt(t1*t1/nu);
  const double sb1 = sqrt(nu/(nu+t1*t1));
  const double a2 = sign(t2)*sqrt(t2*t2/nu);
  const double sb2 = sqrt(nu/(nu+t2*t2));
  size_t j;
  double* C = new double[J];
  for(j=0; j<J; j++){
    double R = sqrt(nu)*(delta1[j] - delta2[j])/(t1-t2);
    double C1 = owent(delta2[j]*sb2, a2) - owent(delta1[j]*sb1, a1);
    double H21 = fabs(t1) < 1 ?
      (t1-(t1-t2)/(1-delta2[j]/delta1[j]))/sqrt(nu) :
      t1/sqrt(nu)*(1-(1-t2/t1)/(1-delta2[j]/delta1[j]));
    double H22 = fabs(t2) < 1 ?
      (t2-(t1-t2)/(delta1[j]/delta2[j]-1))/sqrt(nu) :
      t2/sqrt(nu)*(1-(t1/t2-1)/(delta1[j]/delta2[j]-1));
    double C2 = owent(R, H22) - owent(R, H21);
    //owent(R[j], a2-delta2[j]/R[j) - owent(R[j], a1-delta1[j]/R[j]);
    double H32 = fabs(t2) < 1 ?
      t2/sqrt(nu)*(1- (delta1[j]/delta2[j]-1)/(t1/t2-1)) -
        (delta1[j]/delta2[j]-1)/sqrt(nu)*nu/(t1-t2) :
      t2/sqrt(nu)*(1- (delta1[j]/delta2[j]-1)/(t1/t2-1)) -
        (delta1[j]/delta2[j]-1)/sqrt(nu)*nu/t2/(t1/t2-1);
    double H31 = fabs(t1) < 1 ?
      t1/sqrt(nu)*(1- (1-delta2[j]/delta1[j])/(1-t2/t1)) -
        (1-delta2[j]/delta1[j])/sqrt(nu)*nu/(t1-t2) :
      t1/sqrt(nu)*(1- (1-delta2[j]/delta1[j])/(1-t2/t1)) -
        (1-delta2[j]/delta1[j])/sqrt(nu)*nu/t1/(1-t2/t1);
    double C3 =
      owent(delta2[j]*sb2, H32) - //(delta2[j]*ab2-R[j])/b2/delta2[j]) -
      owent(delta1[j]*sb1, H31);
    C[j] = 2*(C1 - C2 - C3) + (delta1[j] >= 0) - (delta2[j] >= 0);
  }
  return C;
}

// [[Rcpp::export]]
NumericVector OwenCDF4(int nu, double t1, double t2, NumericVector delta1, NumericVector delta2){
  if(t1 <= t2){ // dans R
    return studentCDF(t2, nu, delta2) - studentCDF(t1, nu, delta1);
  }
  const int J = delta1.size();
  NumericVector out(J);
  if(nu == 1){
    double* C = OwenCDF4_C(nu, t1, t2, delta1, delta2, J);
    for(int j=0; j<J; j++){
      out[j] = C[j];
    }
    delete[] C;
    return out;
  }
  const mp::float128 t1t1(t1*t1);
  //const mp::float128 a1 = sign(t1)*mp::sqrt(t1t1/nu);
  const mp::float128 b1 = nu/(nu+t1t1);
  const mp::float128 sb1 = mp::sqrt(b1);
  mp::float128 ab1, asb1;
  if(fabs(t1) > DBL_MAX){ // est-ce utile ?..
    ab1 = mp::float128(0);
    asb1 = mp::float128(sign(t1));
  }else{
    ab1 = mp::float128(sqrt(nu) * 1/(nu/t1+t1));
    asb1 = sign(t1) * mp::sqrt(1/(nu/t1t1+1));
  }
  const mp::float128 t2t2(t2*t2);
  //const mp::float128 a2 = sign(t2)*mp::sqrt(t2t2/nu);
  const mp::float128 b2 = nu/(nu+t2t2);
  const mp::float128 sb2 = mp::sqrt(b2);
  const mp::float128 ab2 = fabs(t2) > DBL_MAX ?
    mp::float128(0) :
    mp::float128(sqrt(nu) * 1/(nu/t2+t2));
  const mp::float128 asb2 = fabs(t2) > DBL_MAX ?
    mp::float128(sign(t2)) :
    sign(t2) * mp::sqrt(1/(nu/t2t2+1));
  mp::float128 R[J];
  mp::float128 dnormdsb1[J];
  mp::float128 dnormdsb2[J];
  mp::float128 Roversb1[J];
  mp::float128 Roversb2[J];
  mp::float128 dabminusRoversb1[J];
  mp::float128 dabminusRoversb2[J];
  mp::float128 dnormR[J];
  mp::float128 RdnormR[J];
  const int n = nu-1;
  mp::float128 M1[n][J];
  mp::float128 M2[n][J];
  mp::float128 H[n][J];
  int j;
  for(j=0; j<J; j++){
    R[j] = mp::float128(sqrt(nu)*(delta1[j]-delta2[j])/(t1-t2));
    dnormdsb1[j] = dnorm128(delta1[j] * sb1);
    dnormdsb2[j] = dnorm128(delta2[j] * sb2);
    Roversb1[j] = fabs(t1) < 1 ?
    R[j]/sb1 :
      sign(t1)*(delta1[j]-delta2[j])*sqrt(nu/t1t1+1)/(1-t2/t1);
    dabminusRoversb1[j] = delta1[j]*asb1 - Roversb1[j];
    Roversb2[j] = fabs(t2) < 1 ?
    R[j]/sb2 :
      sign(t2)*(delta1[j]-delta2[j])*sqrt(nu/t2t2+1)/(t1/t2-1);
    dabminusRoversb2[j] = delta2[j]*asb2 - Roversb2[j];
    dnormR[j] = dnorm128(R[j]);
    RdnormR[j] = R[j] * dnormR[j];
    H[0][j] = -dnormR[j] *
      (pnorm128(asb2*Roversb2[j]-delta2[j]) -
      pnorm128(asb1*Roversb1[j]-delta1[j]));
    M1[0][j] = asb1 * dnormdsb1[j] *
      (pnorm128(delta1[j]*asb1) - pnorm128(dabminusRoversb1[j]));
    M2[0][j] = asb2 * dnormdsb2[j] *
      (pnorm128(delta2[j]*asb2) - pnorm128(dabminusRoversb2[j]));
  }
  if(nu >= 3){
    for(j=0; j<J; j++){
      H[1][j] = R[j] > DBL_MAX ? 0 : R[j] * H[0][j]; // pas besoin car je vais traiter delta1=Inf
      M1[1][j] = delta1[j] > DBL_MAX ? // idem, pas besoin
      0 :
        delta1[j]*ab1*M1[0][j] + ab1 * dnormdsb1[j] *
          (dnorm128(delta1[j]*asb1) - dnorm128(dabminusRoversb1[j]));
      M2[1][j] = delta2[j]*ab2*M2[0][j] + ab2 * dnormdsb2[j] *
        (dnorm128(delta2[j]*asb2) - dnorm128(dabminusRoversb2[j]));
    }
    if(nu >= 4){
      mp::float128 A[n];
      mp::float128 L1[n-2][J];
      mp::float128 L2[n-2][J];
      A[0] = 1;
      A[1] = 1;
      for(j=0; j<J; j++){
        L1[0][j] = ab1 * RdnormR[j] * 0.5*dnorm128(asb1*Roversb1[j]-delta1[j]);
        L2[0][j] = ab2 * RdnormR[j] * 0.5*dnorm128(asb2*Roversb2[j]-delta2[j]);
      }
      int k;
      for(k=2; k<n; k++){
        A[k] = 1.0/k/A[k-1];
      }
      if(nu >= 5){
        for(k=1; k<n-2; k++){
          for(j=0; j<J; j++){
            L1[k][j] = A[k+2] * R[j] * L1[k-1][j];
            L2[k][j] = A[k+2] * R[j] * L2[k-1][j];
          }
        }
      }
      for(k=2; k<n; k++){
        for(j=0; j<J; j++){
          H[k][j] = A[k] * R[j] * H[k-1][j];
          M1[k][j] = (k-1.0)/k *
            (A[k-2] * delta1[j] * ab1 * M1[k-1][j] + b1*M1[k-2][j]) - L1[k-2][j];
          M2[k][j] = (k-1.0)/k *
            (A[k-2] * delta2[j] * ab2 * M2[k-1][j] + b2*M2[k-2][j]) - L2[k-2][j];
        }
      }
    }
  }
  std::vector<mp::float128> sum(J);
  int i;
  if(nu % 2 == 0){
    for(i=0; i<n; i+=2){
      for(j=0; j<J; j++){
        sum[j] += M2[i][j] - M1[i][j] + H[i][j];
      }
    }
    for(j=0; j<J; j++){
      out[j] = pnorm64(-delta2[j]) - pnorm64(-delta1[j]) +
        root_two_pi*sum[j].convert_to<double>();
    }
    return out;
  }else{
    for(i=1; i<n; i+=2){
      for(j=0; j<J; j++){
        sum[j] += M2[i][j] - M1[i][j] + H[i][j];
      }
    }
    double* C = OwenCDF4_C(nu, t1, t2, delta1, delta2, J);
    for(j=0; j<J; j++){
      out[j] = C[j] + 2*sum[j].convert_to<double>();
    }
    delete[] C;
    return out;
  }
}

// --- Owen cumulative function 2 ------------------------------------------- //
double* OwenCDF2_C(int nu, double t1, double t2, NumericVector delta1, NumericVector delta2, size_t J){
  const double sb1 = sqrt(nu/(nu+t1*t1));
  const double sb2 = sqrt(nu/(nu+t2*t2));
  size_t j;
  double* C = new double[J];
  for(j=0; j<J; j++){
    double R = sqrt(nu)*(delta1[j] - delta2[j])/(t1-t2);
    double H21 = fabs(t1) < 1 ?
      (t1-(t1-t2)/(1-delta2[j]/delta1[j]))/sqrt(nu) :
      t1/sqrt(nu)*(1-(1-t2/t1)/(1-delta2[j]/delta1[j]));
    double H22 = fabs(t2) < 1 ?
      (t2-(t1-t2)/(delta1[j]/delta2[j]-1))/sqrt(nu) :
      t2/sqrt(nu)*(1-(t1/t2-1)/(delta1[j]/delta2[j]-1));
    double C2 = owent(R, H22) - owent(R, H21);
    double H32 = fabs(t2) < 1 ?
      t2/sqrt(nu)*(1- (delta1[j]/delta2[j]-1)/(t1/t2-1)) -
        (delta1[j]/delta2[j]-1)/sqrt(nu)*nu/(t1-t2) :
      t2/sqrt(nu)*(1- (delta1[j]/delta2[j]-1)/(t1/t2-1)) -
        (delta1[j]/delta2[j]-1)/sqrt(nu)*nu/t2/(t1/t2-1);
    double H31 = fabs(t1) < 1 ?
      t1/sqrt(nu)*(1- (1-delta2[j]/delta1[j])/(1-t2/t1)) -
        (1-delta2[j]/delta1[j])/sqrt(nu)*nu/(t1-t2) :
      t1/sqrt(nu)*(1- (1-delta2[j]/delta1[j])/(1-t2/t1)) -
        (1-delta2[j]/delta1[j])/sqrt(nu)*nu/t1/(1-t2/t1);
    double C3 =
      owent(delta2[j]*sb2, H32) - owent(delta1[j]*sb1, H31);
    C[j] = -2*(C2 + C3) + (delta1[j] >= 0) - (delta2[j] >= 0) +
      pnorm64(-delta1[j]*sb1) - pnorm64(-delta2[j]*sb2);
  }
  return C;
}

// [[Rcpp::export]]
NumericVector OwenCDF2(int nu, double t1, double t2, NumericVector delta1, NumericVector delta2){
  const int J = delta1.size();
  NumericVector out(J);
  if(nu == 1){
    double* C = OwenCDF2_C(nu, t1, t2, delta1, delta2, J);
    for(int j=0; j<J; j++){
      out[j] = C[j];
    }
    delete[] C;
    return out;
  }
  const mp::float128 t1t1(t1*t1);
  const mp::float128 b1 = nu/(nu+t1t1);
  const mp::float128 sb1 = mp::sqrt(b1);
  const mp::float128 ab1 = mp::float128(sqrt(nu) * 1/(nu/t1+t1));
  const mp::float128 asb1 = sign(t1) * mp::sqrt(1/(nu/t1t1+1));
  const mp::float128 t2t2(t2*t2);
  const mp::float128 b2 = nu/(nu+t2t2);
  const mp::float128 sb2 = mp::sqrt(b2);
  const mp::float128 ab2 = mp::float128(sqrt(nu) * 1/(nu/t2+t2));
  const mp::float128 asb2 = sign(t2) * mp::sqrt(1/(nu/t2t2+1));
  mp::float128 R[J];
  mp::float128 dnormdsb1[J];
  mp::float128 dnormdsb2[J];
  mp::float128 Roversb1[J];
  mp::float128 Roversb2[J];
  mp::float128 dabminusRoversb1[J];
  mp::float128 dabminusRoversb2[J];
  mp::float128 dnormR[J];
  mp::float128 RdnormR[J];
  const int n = nu-1;
  mp::float128 M1[n][J];
  mp::float128 M2[n][J];
  mp::float128 H[n][J];
  int j;
  for(j=0; j<J; j++){
    R[j] = mp::float128(sqrt(nu)*(delta1[j]-delta2[j])/(t1-t2));
    dnormdsb1[j] = dnorm128(delta1[j] * sb1);
    dnormdsb2[j] = dnorm128(delta2[j] * sb2);
    Roversb1[j] = fabs(t1) < 1 ?
      R[j]/sb1 :
      sign(t1)*(delta1[j]-delta2[j])*sqrt(nu/t1t1+1)/(1-t2/t1);
    dabminusRoversb1[j] = delta1[j]*asb1 - Roversb1[j];
    Roversb2[j] = fabs(t2) < 1 ?
      R[j]/sb2 :
      sign(t2)*(delta1[j]-delta2[j])*sqrt(nu/t2t2+1)/(t1/t2-1);
    dabminusRoversb2[j] = delta2[j]*asb2 - Roversb2[j];
    dnormR[j] = dnorm128(R[j]);
    RdnormR[j] = R[j] * dnormR[j];
    H[0][j] = -dnormR[j] *
      (pnorm128(asb2*Roversb2[j]-delta2[j]) -
        pnorm128(asb1*Roversb1[j]-delta1[j]));
    M1[0][j] = asb1 * dnormdsb1[j] * pnorm128(dabminusRoversb1[j]);
    M2[0][j] = asb2 * dnormdsb2[j] * pnorm128(dabminusRoversb2[j]);
  }
  if(nu >= 3){
    for(j=0; j<J; j++){
      H[1][j] = R[j] * H[0][j];
      M1[1][j] = delta1[j]*ab1*M1[0][j] + ab1 * dnormdsb1[j] *
                  dnorm128(dabminusRoversb1[j]);
      M2[1][j] = delta2[j]*ab2*M2[0][j] + ab2 * dnormdsb2[j] *
                  dnorm128(dabminusRoversb2[j]);
    }
    if(nu >= 4){
      mp::float128 A[n];
      mp::float128 L1[n-2][J];
      mp::float128 L2[n-2][J];
      A[0] = 1;
      A[1] = 1;
      for(j=0; j<J; j++){
        L1[0][j] = ab1 * RdnormR[j] * 0.5*dnorm128(asb1*Roversb1[j]-delta1[j]);
        L2[0][j] = ab2 * RdnormR[j] * 0.5*dnorm128(asb2*Roversb2[j]-delta2[j]);
      }
      int k;
      for(k=2; k<n; k++){
        A[k] = 1.0/k/A[k-1];
      }
      if(nu >= 5){
        for(k=1; k<n-2; k++){
          for(j=0; j<J; j++){
            L1[k][j] = A[k+2] * R[j] * L1[k-1][j];
            L2[k][j] = A[k+2] * R[j] * L2[k-1][j];
          }
        }
      }
      for(k=2; k<n; k++){
        for(j=0; j<J; j++){
          H[k][j] = A[k] * R[j] * H[k-1][j];
          M1[k][j] = (k-1.0)/k *
            (A[k-2] * delta1[j] * ab1 * M1[k-1][j] + b1*M1[k-2][j]) + L1[k-2][j];
          M2[k][j] = (k-1.0)/k *
            (A[k-2] * delta2[j] * ab2 * M2[k-1][j] + b2*M2[k-2][j]) + L2[k-2][j];
        }
      }
    }
  }
  std::vector<mp::float128> sum(J);
  int i;
  if(nu % 2 == 0){
    for(i=0; i<n; i+=2){
      for(j=0; j<J; j++){
        sum[j] += M1[i][j] - M2[i][j] - H[i][j];
      }
    }
    for(j=0; j<J; j++){
      out[j] = root_two_pi*sum[j].convert_to<double>();
    }
    return out;
  }else{
    for(i=1; i<n; i+=2){
      for(j=0; j<J; j++){
        sum[j] += M1[i][j] - M2[i][j] - H[i][j];
      }
    }
    double* C = OwenCDF2_C(nu, t1, t2, delta1, delta2, J);
    for(j=0; j<J; j++){
      out[j] = C[j] + 2*sum[j].convert_to<double>();
    }
    delete[] C;
    return out;
  }
}
