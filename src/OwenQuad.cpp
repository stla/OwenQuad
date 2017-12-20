// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
using namespace Rcpp;
#include <boost/multiprecision/float128.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/owens_t.hpp>
#include <cmath>
//#include <climits>

namespace mp = boost::multiprecision;
namespace m = boost::math;

// double AA(int k){
//   if(k ==0){
//     return 1.0;
//   }
//   if(k % 2 == 1){
//     return pow(2.0, k-1) * pow(tgamma((k+1)/2), 2) / tgamma(k+1);
//   }else{
//     return tgamma(k) / pow(2.0, k-2) / pow(tgamma(k/2), 2) / k;
//   }
// }

// double logAA(int k){
//   if(k == 0){
//     return 0.0;
//   }
//   if(k % 2 == 1){
//     return (k-1)*log(2.0) + 2* lgamma((k+1)/2.0) - lgamma(k+1);
//   }else{
//     return lgamma(k) - (k-2)*log(2.0) - 2.0*lgamma(k/2.0) - log(k);
//   }
// }

// mp::float128 logAA(int k){
//   if(k == 0){
//     return 0.0;
//   }
//   mp::float128 kk = mp::float128(k);
//   if(k % 2 == 1){
//     return (k-1)*mp::log(mp::float128(2.0)) + 2*mp::lgamma((kk+1)/2.0) - mp::lgamma(kk+1);
//   }else{
//     return mp::lgamma(kk) - (k-2)*mp::log(mp::float128(2.0)) - 2.0*mp::lgamma(kk/2.0) - mp::log(kk);
//   }
// }

// mp::float128 AA(int k){
//   return mp::exp(logAA(k));
// }
 


// double test(int k){
//   mp::float128 x = ldoublefact(k);
//   return x.convert_to<double>();
// }

const double one_div_root_two = m::constants::one_div_root_two<double>();
const double root_two_pi = m::constants::root_two_pi<double>();
const mp::float128 root_two_pi128 = m::constants::root_two_pi<mp::float128>();
const mp::float128 one_div_root_two_pi128 = m::constants::one_div_root_two_pi<mp::float128>();
const mp::float128 one_div_root_two128 = m::constants::one_div_root_two<mp::float128>();
const mp::float128 log_root_two_pi128 = m::constants::log_root_two_pi<mp::float128>();
const mp::float128 log_two128 = 
  0.69314718055994530941723212145817656807550013436025525412068000949339362196969471560586332699641868754200148102Q;

double pnorm64(double q){
  if(std::isnan(q)){
    return nan("");
  }
  return m::erfc(-q * one_div_root_two)/2.0;
}

mp::float128 dnorm128(mp::float128 x){
  return mp::exp(-x*x/2) * one_div_root_two_pi128;
}

mp::float128 xdnormx(mp::float128 x){
  return mp::exp(mp::log(x) - 0.5*x*x - log_root_two_pi128);
}

mp::float128 xndnorm(mp::float128 x, int n){
  return mp::exp(n*mp::log(x) - 0.5*x*x - log_root_two_pi128);
}

mp::float128 ldoublefact(int k){
  const mp::float128 kk = mp::float128(k);
  if(k % 2 == 0){
    return 0.5*kk*log_two128 + m::lgamma(1.0+0.5*kk);
  }else{
    return m::lgamma(kk+1.0) - 0.5*(kk-1.0)*log_two128 - m::lgamma((kk+1.0)/2.0);
  }
}
// double ldoublefact(int k){
//   if(k % 2 == 0){
//     return k/2*log(2.0) + lgamma(1.0+k/2.0);
//   }else{
//     return lgamma(k+1.0) - ((k-1.0)/2.0)*log(2.0) - lgamma((k+1.0)/2.0);
//   }
// }

mp::float128 anxndnorm(mp::float128 x, int n, int i=0){
  return mp::exp(-ldoublefact(n+i) + n*mp::log(x) - 0.5*x*x - log_root_two_pi128);
}

mp::float128 pnorm128(mp::float128 q){
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
NumericVector studentCDF_C(double q, int nu, NumericVector delta){
  const double a = sign(q)*sqrt(q*q/nu);
  const double sb = sqrt(nu/(nu+q*q));
  const size_t J = delta.size();
  NumericVector C(J);
  for(size_t j=0; j<J; j++){
    C[j] = 2*owent(delta[j] * sb, a) + pnorm64(-delta[j]*sb);
  }
  return C;
}

// [[Rcpp::export]]
NumericVector studentCDF(double q, size_t nu, NumericVector delta){
  if(nu==1){
    return studentCDF_C(q, nu, delta);
  }
  const size_t J = delta.size();
  size_t j;
  NumericVector out(J);
  const mp::float128 qq(q*q);
  const mp::float128 a = sign(q)*mp::sqrt(qq/nu);
  const mp::float128 b = nu/(nu+qq);
  const mp::float128 sb = mp::sqrt(b);
  mp::float128 M[nu-1][J];
  mp::float128 dsb;
  for(j=0; j<J ; j++){
    dsb = delta[j] * sb;
    M[0][j] = a * sb * dnorm128(dsb) * pnorm128(a*dsb);
  }
  if(nu>2){
    for(j=0; j<J; j++){
      M[1][j] = b * (delta[j] * a * M[0][j] + a * dnorm128(delta[j]) * one_div_root_two_pi128);
    }
    if(nu>3){
      mp::float128 A[nu-3]; A[0] = 1.0;
      size_t k;
      if(nu>4){
        for(k=1; k<nu-3; k++){
          A[k] = 1.0/(k*A[k-1]);
        }
      }
      for(k=2; k<nu-1; k++){
        for(j=0; j<J; j++){
          M[k][j] = (k-1) * b * (A[k-2] * delta[j] * a * M[k-1][j] + M[k-2][j]) / k;
        }
      }
    }
  }
  if(nu%2==1){
    NumericVector C = studentCDF_C(q, nu, delta);
    for(j=0; j<J; j++){
      mp::float128 sum = 0;
      for(size_t i=1; i<nu-1; i+=2){
        sum += M[i][j];
      }
      out[j] = C[j] + 2*sum.convert_to<double>();
    }
  }else{
    for(j=0; j<J; j++){
      mp::float128 sum=0; mp::float128 out128;
      for(size_t i=0; i<nu-1; i+=2){
        sum += M[i][j];
      }
      out128 = pnorm128(-delta[j]) + root_two_pi128*sum;
      out[j] = out128.convert_to<double>();
    }
  }
  return out;
}


// ------- Owen first Q-function -------------------------------------------- //
NumericVector OwenQ1_C(int nu, double t, NumericVector delta, NumericVector R){
  const double a = sign(t)*sqrt(t*t/nu);
  const double b = nu/(nu+t*t);
  const double sb = sqrt(b);
  const double ab = sqrt(nu)/(nu/t+t);
  size_t J = delta.size();
  NumericVector C(J);
  for(size_t i=0; i<J; i++){
    double C1 = owent(delta[i]*sb, a);
    double C2 = owent(R[i], a-delta[i]/R[i]);
    double C3 = owent(delta[i]*sb, (ab-R[i]/delta[i])/b);
    C[i] = pnorm64(R[i]) - (delta[i] >= 0) + 2*(C1 - C2 - C3);
  }
  return C;
}

// [[Rcpp::export]]
NumericVector OwenQ1(size_t nu, double t, NumericVector delta, NumericVector R, int algo=1){
  if(nu == 1){
    return OwenQ1_C(nu, t, delta, R);
  }
  const size_t J = delta.size();
  size_t j;
  NumericVector out(J);
  const mp::float128 tt(t*t);
  const mp::float128 a = sign(t)*mp::sqrt(tt/nu);
  const mp::float128 b = nu/(nu+tt);
  const mp::float128 sb = mp::sqrt(b);
  const mp::float128 ab = mp::float128(sqrt(nu)/(nu/t+t));
  const mp::float128 asb = sign(t)/mp::sqrt(nu/tt+1);
  mp::float128 dnormdsb[J];
  mp::float128 dabminusRoversb[J];
  const size_t n = nu-1;
  mp::float128 H[n][J]; mp::float128 M[n][J];
  mp::float128 Lfactor[J];
  for(j=0; j<J; j++){
    dnormdsb[j] = dnorm128(delta[j] * sb);
    dabminusRoversb[j] = (delta[j]*ab - R[j])/sb;
    Lfactor[j] = ab * dnorm128(a*R[j]-delta[j]);
    H[0][j] = dnorm128(R[j]);
    M[0][j] = asb * dnormdsb[j] * (pnorm128(delta[j]*asb) - pnorm128(dabminusRoversb[j]));
  }
  if(nu >= 3){
    for(j=0; j<J; j++){
      H[1][j] = xdnormx(R[j]);
      M[1][j] = delta[j]*ab*M[0][j] + ab * dnormdsb[j] *
                  (dnorm128(delta[j]*asb) - dnorm128(dabminusRoversb[j]));
    }
    if(nu >= 4){
      size_t k;
      if(algo == 1){
        mp::float128 A[n]; A[0] = 1; A[1] = 1;
        mp::float128 L[n-2][J];
        for(j=0; j<J; j++){
          L[0][j] = 0.5*H[1][j];
        }
        for(k=2; k<n; k++){
          A[k] = 1.0/(k*A[k-1]);
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
            M[k][j] = (k-1.0)/k * (A[k-2] * delta[j] * ab * M[k-1][j] + b*M[k-2][j]) - 
                        Lfactor[j]*L[k-2][j];
          }
        }
      }else{ // algo 2
        mp::float128 A[n-1]; A[0] = mp::float128(1);
        mp::float128 halfRR[J]; mp::float128 logR[j];
        for(j=0; j<J; j++){
          mp::float128 Rj = mp::float128(R[j]);
          halfRR[j] = 0.5*Rj*Rj;
          logR[j] = mp::log(Rj);
        }
        for(k=0; k<n-2; k++){
          A[k+1] = 1.0/((k+1)*A[k]); // un de trop
          mp::float128 ldf;
          if(k % 2 == 0){
            ldf = 0.5*(k+2)*log_two128 + m::lgamma(mp::float128(0.5*k+2.0));
          }else{
            ldf = m::lgamma(mp::float128(k+3)) - 0.5*(mp::float128(k+1))*log_two128 - 
                    m::lgamma(mp::float128(0.5*(k+3)));
          }
          mp::float128 r = mp::float128(k+1)/mp::float128(k+2);
          for(j=0; j<J; j++){
            mp::float128 K =
              mp::exp(-ldf + (k+1)*logR[j] - halfRR[j] - log_root_two_pi128);
            H[k+2][j] = K*R[j];
            M[k+2][j] = r *
              (A[k] * delta[j] * ab * M[k+1][j] + b*M[k][j]) - K*Lfactor[j];
          }
        }
      }
    }
  }
  if(nu % 2 == 0){
    for(j=0; j<J; j++){
      mp::float128 sumH = 0; mp::float128 sumM = 0;
      for(size_t i=0; i<n; i+=2){
        sumH += H[i][j];
        sumM += M[i][j];
      }
      mp::float128 out128 = root_two_pi128*(sumM - pnorm128(a*R[j]-delta[j])*sumH) + 
                pnorm128(-delta[j]); 
      out[j] = out128.convert_to<double>();
    }
    return out;
  }else{
    NumericVector C = OwenQ1_C(nu, t, delta, R);
    for(j=0; j<J; j++){
      mp::float128 sumH = 0; mp::float128 sumM = 0;
      for(size_t i=1; i<n; i+=2){
        sumH += H[i][j];
        sumM += M[i][j];
      }
      mp::float128 out128 = 2*(sumM - pnorm128(a*R[j]-delta[j])*sumH);
      out[j] = C[j] + out128.convert_to<double>();
    }
    return out;
  }
}

// --- Owen second Q-function ------------------------------------------- //
NumericVector OwenQ2_C(int nu, double t, NumericVector delta, NumericVector R){
  const double a = sign(t)*sqrt(t*t/nu);
  const double b = nu/(nu+t*t);
  const double ab = sqrt(nu) * 1/(nu/t+t);
  const double sb = sqrt(nu/(nu+t*t));
  size_t J = delta.size();
  NumericVector C(J);
  for(size_t j=0; j<J; j++){
    double C2 = owent(R[j], a-delta[j]/R[j]);
    double C3 = owent(delta[j]*sb, (ab-R[j]/delta[j])/b);
    C[j] = 2*(C2 + C3) + (delta[j] >= 0) +
            pnorm64(-delta[j]*sb) - pnorm64(R[j]);
  }
  return C;
}

// [[Rcpp::export]]
NumericVector OwenQ2
    (size_t nu, double t, NumericVector delta, NumericVector R, int algo=1){
  if(nu == 1){
    return OwenQ2_C(nu, t, delta, R);
  }
  const size_t J = delta.size();
  NumericVector out(J);
  const mp::float128 tt(t*t);
  const mp::float128 a = sign(t)*mp::sqrt(tt/nu);
  const mp::float128 b = nu/(nu+tt);
  const mp::float128 sb = mp::sqrt(b);
  const mp::float128 ab = mp::float128(sqrt(nu)/(nu/t+t));
  const mp::float128 asb = sign(t)/mp::sqrt(nu/tt+1);
  mp::float128 dnormdsb[J];
  mp::float128 dabminusRoversb[J];
  const size_t n = nu-1;
  mp::float128 H[n][J]; mp::float128 M[n][J];
  mp::float128 Lfactor[J];
  size_t j;
  for(j=0; j<J; j++){
    dnormdsb[j] = dnorm128(delta[j] * sb);
    dabminusRoversb[j] = (delta[j]*ab - R[j])/sb;
    H[0][j] = dnorm128(R[j]);
    M[0][j] = asb * dnormdsb[j] * pnorm128(dabminusRoversb[j]);
    Lfactor[j] = ab * dnorm128(a*R[j]-delta[j]);
  }
  if(nu >= 3){
    for(j=0; j<J; j++){
      H[1][j] = xdnormx(R[j]);
      M[1][j] = delta[j]*ab*M[0][j] + ab * dnormdsb[j] *
                  dnorm128(dabminusRoversb[j]);
    }
    if(nu >= 4){
      size_t k;
      if(algo == 1){
        mp::float128 A[n]; A[0] = 1; A[1] = 1;
        mp::float128 L[n-2][J];
        for(j=0; j<J; j++){
          L[0][j] = 0.5*H[1][j];
        }
        for(k=2; k<n; k++){
          A[k] = 1.0/(k*A[k-1]);
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
            M[k][j] = (k-1.0)/k * (A[k-2] * delta[j] * ab * M[k-1][j] + b*M[k-2][j]) + 
                        Lfactor[j]*L[k-2][j];
          }
        }
      }else{ // algo 2
        mp::float128 A[n-1]; A[0] = mp::float128(1);
        mp::float128 halfRR[J]; mp::float128 logR[j];
        for(j=0; j<J; j++){
          mp::float128 Rj = mp::float128(R[j]);
          halfRR[j] = 0.5*Rj*Rj;
          logR[j] = mp::log(Rj);
        }
        for(k=0; k<n-2; k++){
          A[k+1] = 1.0/((k+1)*A[k]); // un de trop
          mp::float128 ldf;
          if(k % 2 == 0){
            ldf = 0.5*(k+2)*log_two128 + m::lgamma(mp::float128(0.5*k+2.0));
          }else{
            ldf = m::lgamma(mp::float128(k+3)) - 0.5*(mp::float128(k+1))*log_two128 - 
                    m::lgamma(mp::float128(0.5*(k+3)));
          }
          mp::float128 r = mp::float128(k+1)/mp::float128(k+2);
          for(j=0; j<J; j++){
            mp::float128 K =
              mp::exp(-ldf + (k+1)*logR[j] - halfRR[j] - log_root_two_pi128);
            H[k+2][j] = K*R[j];
            M[k+2][j] = r *
              (A[k] * delta[j] * ab * M[k+1][j] + b*M[k][j]) + K*Lfactor[j];
          } 
        }
      }
    }
  }
  if(nu % 2 == 0){
    for(j=0; j<J; j++){
      mp::float128 sumH = 0; mp::float128 sumM = 0;
      for(size_t i=0; i<n; i+=2){
        sumH += H[i][j];
        sumM += M[i][j];
      }
      mp::float128 out128 = root_two_pi128*(sumM + pnorm128(a*R[j]-delta[j])*sumH); 
      out[j] = out128.convert_to<double>();
    }
    return out;
  }else{
    NumericVector C = OwenQ2_C(nu, t, delta, R);
    for(j=0; j<J; j++){
      mp::float128 sumH = 0; mp::float128 sumM = 0;
      for(size_t i=1; i<n; i+=2){
        sumH += H[i][j];
        sumM += M[i][j];
      }
      mp::float128 out128 = 2*(sumM + pnorm128(a*R[j]-delta[j])*sumH);
      out[j] = C[j] + out128.convert_to<double>();
    }
    return out;
  }
}

// --- Owen cumulative function 4 ------------------------------------------- //
NumericVector OwenCDF4_C(int nu, double t1, double t2,
                                 NumericVector delta1, NumericVector delta2){
  const double a1 = sign(t1)*sqrt(t1*t1/nu);
  const double sb1 = sqrt(nu/(nu+t1*t1));
  const double a2 = sign(t2)*sqrt(t2*t2/nu);
  const double sb2 = sqrt(nu/(nu+t2*t2));
  size_t J = delta1.size();
  size_t j;
  NumericVector C(J);
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
NumericVector OwenCDF4(size_t nu, double t1, double t2, NumericVector delta1, 
                                              NumericVector delta2, int algo=1){
  if(nu == 1){
    return OwenCDF4_C(nu, t1, t2, delta1, delta2);
  }
  const size_t J = delta1.size();
  size_t j;
  NumericVector out(J);
  const mp::float128 t1t1(t1*t1);
  const mp::float128 b1 = nu/(nu+t1t1);
  const mp::float128 sb1 = mp::sqrt(b1);
  const mp::float128 ab1 = mp::float128(sqrt(nu)/(nu/t1+t1));
  const mp::float128 asb1 = sign(t1)/mp::sqrt(nu/t1t1+1);
  const mp::float128 t2t2(t2*t2);
  const mp::float128 b2 = nu/(nu+t2t2);
  const mp::float128 sb2 = mp::sqrt(b2);
  const mp::float128 ab2 = mp::float128(sqrt(nu)/(nu/t2+t2));
  const mp::float128 asb2 = sign(t2)/mp::sqrt(nu/t2t2+1);
  mp::float128 R[J];
  mp::float128 dnormdsb1[J]; mp::float128 dnormdsb2[J];
  mp::float128 aRminusdelta1[J]; mp::float128 aRminusdelta2[J];
  mp::float128 dabminusRoversb1[J]; mp::float128 dabminusRoversb2[J];
  const size_t n = nu-1;
  mp::float128 M1[n][J]; mp::float128 M2[n][J]; mp::float128 H[n][J];
  mp::float128 Lfactor1[J]; mp::float128 Lfactor2[J];
  for(j=0; j<J; j++){
    R[j] = mp::float128(sqrt(nu)*(delta1[j]-delta2[j])/(t1-t2));
    dnormdsb1[j] = dnorm128(delta1[j] * sb1);
    dnormdsb2[j] = dnorm128(delta2[j] * sb2);
    mp::float128 Roversb1 = fabs(t1) < 1 ?
      R[j]/sb1 :
      sign(t1)*(delta1[j]-delta2[j])*mp::sqrt(nu/t1t1+1)/(1-t2/t1);
    dabminusRoversb1[j] = delta1[j]*asb1 - Roversb1;
    mp::float128 Roversb2 = fabs(t2) < 1 ?
      R[j]/sb2 :
      sign(t2)*(delta1[j]-delta2[j])*mp::sqrt(nu/t2t2+1)/(t1/t2-1);
    dabminusRoversb2[j] = delta2[j]*asb2 - Roversb2;
    H[0][j] = dnorm128(R[j]);
    M1[0][j] = asb1 * dnormdsb1[j] *
                (pnorm128(delta1[j]*asb1) - pnorm128(dabminusRoversb1[j]));
    M2[0][j] = asb2 * dnormdsb2[j] *
                (pnorm128(delta2[j]*asb2) - pnorm128(dabminusRoversb2[j]));
    aRminusdelta1[j] = asb1*Roversb1-delta1[j];
    aRminusdelta2[j] = asb2*Roversb2-delta2[j];
    Lfactor1[j] = ab1 * dnorm128(aRminusdelta1[j]);
    Lfactor2[j] = ab2 * dnorm128(aRminusdelta2[j]);
  }
  if(nu >= 3){
    for(j=0; j<J; j++){
      H[1][j] = xdnormx(R[j]);
      M1[1][j] = delta1[j]*ab1*M1[0][j] + ab1 * dnormdsb1[j] *
                  (dnorm128(delta1[j]*asb1) - dnorm128(dabminusRoversb1[j]));
      M2[1][j] = delta2[j]*ab2*M2[0][j] + ab2 * dnormdsb2[j] *
                  (dnorm128(delta2[j]*asb2) - dnorm128(dabminusRoversb2[j]));
    }
    if(nu >= 4){
      size_t k;
      if(algo == 1){
        mp::float128 A[n]; A[0] = 1; A[1] = 1;
        mp::float128 L[n-2][J];
        for(j=0; j<J; j++){
          L[0][j] = 0.5*H[1][j]; 
        }
        for(k=2; k<n; k++){
          A[k] = 1.0/(k*A[k-1]);
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
            M1[k][j] = (k-1.0)/k *
              (A[k-2] * delta1[j] * ab1 * M1[k-1][j] + b1*M1[k-2][j]) - 
                Lfactor1[j]*L[k-2][j];
            M2[k][j] = (k-1.0)/k *
              (A[k-2] * delta2[j] * ab2 * M2[k-1][j] + b2*M2[k-2][j]) - 
                Lfactor2[j]*L[k-2][j];
          }
        }
      }else{ // algo2
        mp::float128 A[n-1]; A[0] = mp::float128(1);
        mp::float128 halfRR[J]; mp::float128 logR[j];
        for(j=0; j<J; j++){
          mp::float128 Rj = mp::float128(R[j]);
          halfRR[j] = 0.5*Rj*Rj;
          logR[j] = mp::log(Rj);
        }
        for(k=0; k<n-2; k++){
          A[k+1] = 1.0/((k+1)*A[k]); // un de trop
          mp::float128 ldf;
          if(k % 2 == 0){
            ldf = 0.5*(k+2)*log_two128 + m::lgamma(mp::float128(0.5*k+2.0));
          }else{
            ldf = m::lgamma(mp::float128(k+3)) - 0.5*(mp::float128(k+1))*log_two128 - 
                    m::lgamma(mp::float128(0.5*(k+3)));
          }
          mp::float128 r = mp::float128(k+1)/mp::float128(k+2);
          for(j=0; j<J; j++){
            mp::float128 K =
              mp::exp(-ldf + (k+1)*logR[j] - halfRR[j] - log_root_two_pi128);
            H[k+2][j] = K*R[j];
            M1[k+2][j] = r *
              (A[k] * delta1[j] * ab1 * M1[k+1][j] + b1*M1[k][j]) -
                K*Lfactor1[j];
            M2[k+2][j] = r *
              (A[k] * delta2[j] * ab2 * M2[k+1][j] + b2*M2[k][j]) -
                K*Lfactor2[j];
          }
        }
      }
    }
  }
  if(nu % 2 == 0){
    for(j=0; j<J; j++){
      mp::float128 sumH=0; mp::float128 sumM=0; mp::float128 out128;
      for(size_t i=0; i<n; i+=2){
        sumH += H[i][j];
        sumM += M2[i][j] - M1[i][j];
      }
      out128 = root_two_pi128*(sumM + 
                (pnorm128(aRminusdelta1[j]) - pnorm128(aRminusdelta2[j]))*sumH) + 
                pnorm128(-delta2[j]) - pnorm128(-delta1[j]);
      out[j] = out128.convert_to<double>();
    }
    return out;
  }else{
    NumericVector C = OwenCDF4_C(nu, t1, t2, delta1, delta2);
    for(j=0; j<J; j++){
      mp::float128 sumH=0; mp::float128 sumM=0; mp::float128 out128;
      for(size_t i=1; i<n; i+=2){
        sumH += H[i][j];
        sumM += M2[i][j] - M1[i][j];
      }
      out128 = 2*(sumM + 
                (pnorm128(aRminusdelta1[j]) - pnorm128(aRminusdelta2[j]))*sumH);
      out[j] = C[j] + out128.convert_to<double>();
    }
    return out;
  }
}

// NumericVector oo(size_t nu, double t1, double t2, NumericVector delta1, NumericVector delta2){
//   const size_t J = delta1.size();
//   size_t j;
//   NumericVector out(J);
//   const mp::float128 t1t1(t1*t1);
//   const mp::float128 nuqp = mp::float128(nu);
//   const mp::float128 b1 = nuqp/(nuqp+t1t1);
//   const mp::float128 sb1 = mp::sqrt(b1);
//   const mp::float128 ab1 = mp::sqrt(nuqp)/(nuqp/t1+t1);
//   const mp::float128 asb1 = sign(t1) * mp::sqrt(1/(nuqp/t1t1+1));
//   const mp::float128 t2t2(t2*t2);
//   const mp::float128 b2 = nuqp/(nuqp+t2t2);
//   const mp::float128 sb2 = mp::sqrt(b2);
//   const mp::float128 ab2 = mp::sqrt(nuqp)/(nuqp/t2+t2);
//   const mp::float128 asb2 = sign(t2) * mp::sqrt(1/(nuqp/t2t2+1));
//   mp::float128 R[J];
//   mp::float128 dnormdsb1[J];
//   mp::float128 dnormdsb2[J];
//   mp::float128 Roversb1[J];
//   mp::float128 Roversb2[J];
//   mp::float128 dabminusRoversb1[J];
//   mp::float128 dabminusRoversb2[J];
//   const size_t n = nu-1;
//   mp::float128 M1[n][J];
//   mp::float128 M2[n][J];
//   mp::float128 H[n][J];
//   mp::float128 Hfactor[J];
//   mp::float128 G1[J];
//   mp::float128 G2[J];
//   for(j=0; j<J; j++){
//     R[j] = mp::sqrt(nuqp)*(delta1[j]-delta2[j])/(t1-t2);
//     dnormdsb1[j] = dnorm128(delta1[j] * sb1);
//     dnormdsb2[j] = dnorm128(delta2[j] * sb2);
//     Roversb1[j] = fabs(t1) < 1 ?
//       R[j]/sb1 :
//       sign(t1)*(delta1[j]-delta2[j])*mp::sqrt(nuqp/t1t1+1)/(1-t2/t1);
//     dabminusRoversb1[j] = delta1[j]*asb1 - Roversb1[j];
//     Roversb2[j] = fabs(t2) < 1 ?
//       R[j]/sb2 :
//       sign(t2)*(delta1[j]-delta2[j])*mp::sqrt(nu/t2t2+1)/(t1/t2-1);
//     dabminusRoversb2[j] = delta2[j]*asb2 - Roversb2[j];
//     mp::float128 c1 = asb1*Roversb1[j]-delta1[j];
//     mp::float128 c2 = asb2*Roversb2[j]-delta2[j];
//     Hfactor[j] = pnorm128(c2) - pnorm128(c1);
//     H[0][j] = -dnorm128(R[j]);
//     M1[0][j] = asb1 * dnormdsb1[j] *
//       (pnorm128(delta1[j]*asb1) - pnorm128(dabminusRoversb1[j]));
//     M2[0][j] = asb2 * dnormdsb2[j] *
//       (pnorm128(delta2[j]*asb2) - pnorm128(dabminusRoversb2[j]));
//     G1[j] = ab1*dnorm128(c1); // utilisé seulement si nu >= 5
//     G2[j] = ab2*dnorm128(c2); // ''
//   }
//   if(nu >= 3){
//     for(j=0; j<J; j++){
//       H[1][j] = -xndnorm(R[j],1);
//       M1[1][j] = delta1[j]*ab1*M1[0][j] + ab1 * dnormdsb1[j] *
//         (dnorm128(delta1[j]*asb1) - dnorm128(dabminusRoversb1[j]));
//       M2[1][j] = delta2[j]*ab2*M2[0][j] + ab2 * dnormdsb2[j] *
//         (dnorm128(delta2[j]*asb2) - dnorm128(dabminusRoversb2[j]));
//     }
//     if(nu >= 4){
//       mp::float128 A[n-1];
//       A[0] = mp::float128(1);
//       size_t k;
//       for(k=0; k<n-2; k++){
//         A[k+1] = 1.0/((k+1)*A[k]); // un de trop
//         mp::float128 k2 = mp::float128(k+2);
//         mp::float128 ldf;
//         if(k % 2 == 0){
//           ldf = 0.5*k2*log_two128 + m::lgamma(1.0+0.5*k2);
//         }else{
//           ldf = m::lgamma(k2+1.0) - 0.5*(k2-1.0)*log_two128 - m::lgamma((k2+1.0)/2.0);
//         }
//         mp::float128 r = mp::float128(k+1)/k2;
//         for(j=0; j<J; j++){
//           //mp::float128 K = anxndnorm(R[j],k+1,1); // sortir ce qui est en k
//           mp::float128 K =
//             mp::exp(-ldf + (k+1)*mp::log(R[j]) - 0.5*R[j]*R[j] - log_root_two_pi128);
//           H[k+2][j] = -R[j]*K;// = -anxndnorm(R[j],k+2);
//           M1[k+2][j] = r *
//             (A[k] * delta1[j] * ab1 * M1[k+1][j] + b1*M1[k][j]) -
//               K*G1[j];
//           M2[k+2][j] = r *
//             (A[k] * delta2[j] * ab2 * M2[k+1][j] + b2*M2[k][j]) -
//               K*G2[j];
//         }
//       }
//     }
//   }
//   // NumericVector m1(n);
//   // NumericVector m2(n);
//   // NumericVector h(n);
//   // for(int k=0; k<n; k++){
//   //   m1[k] = M1[k][1].convert_to<double>();
//   //   m2[k] = M2[k][1].convert_to<double>();
//   //   h[k] = H[k][1].convert_to<double>();
//   // }
//   // List ret; ret["m1"] = m1; ret["m2"] = m2; ret["h"] = h;
//   // return ret;
//   size_t i;
//   if(nu % 2 == 0){
//     for(j=0; j<J; j++){
//       mp::float128 sumH, sumM, sum;
//       for(i=0; i<n; i+=2){
//         sumH += H[i][j];
//         sumM += M2[i][j] - M1[i][j];
//       }
//       sum = sumM + Hfactor[j]*sumH; // faire out128 ?
//       out[j] = pnorm64(-delta2[j]) - pnorm64(-delta1[j]) +
//         root_two_pi*sum.convert_to<double>();
//     }
//     return out;
//   }else{
//     std::vector<mp::float128> sum(J);
//     for(i=1; i<n; i+=2){
//       for(j=0; j<J; j++){
//         sum[j] += M2[i][j] - M1[i][j] + Hfactor[j]*H[i][j];
//       }
//     }
//     NumericVector C = OwenCDF4_C(nu, t1, t2, delta1, delta2);
//     for(j=0; j<J; j++){
//       out[j] = C[j] + 2*sum[j].convert_to<double>();
//     }
//     return out;
//   }
// }


// --- Owen cumulative function 2 ------------------------------------------- //
NumericVector OwenCDF2_C
    (int nu, double t1, double t2, NumericVector delta1, NumericVector delta2){
  const double sb1 = sqrt(nu/(nu+t1*t1));
  const double sb2 = sqrt(nu/(nu+t2*t2));
  const size_t J = delta1.size();
  NumericVector C(J);
  for(size_t j=0; j<J; j++){
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
NumericVector OwenCDF2
    (size_t nu, double t1, double t2, NumericVector delta1, NumericVector delta2, int algo=1){
  if(nu == 1){
    return OwenCDF2_C(nu, t1, t2, delta1, delta2);
  }
  const size_t J = delta1.size();
  size_t j;
  NumericVector out(J);
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
  mp::float128 dnormdsb1[J]; mp::float128 dnormdsb2[J];
  mp::float128 aRminusdelta1[J]; mp::float128 aRminusdelta2[J];
  mp::float128 dabminusRoversb1[J]; mp::float128 dabminusRoversb2[J];
  const size_t n = nu-1;
  mp::float128 M1[n][J]; mp::float128 M2[n][J]; mp::float128 H[n][J];
  mp::float128 Lfactor1[J]; mp::float128 Lfactor2[J];  
  for(j=0; j<J; j++){
    R[j] = mp::float128(sqrt(nu)*(delta1[j]-delta2[j])/(t1-t2));
    dnormdsb1[j] = dnorm128(delta1[j] * sb1);
    dnormdsb2[j] = dnorm128(delta2[j] * sb2);
    mp::float128 Roversb1 = fabs(t1) < 1 ?
      R[j]/sb1 :
      sign(t1)*(delta1[j]-delta2[j])*mp::sqrt(nu/t1t1+1)/(1-t2/t1);
    dabminusRoversb1[j] = delta1[j]*asb1 - Roversb1;
    mp::float128 Roversb2 = fabs(t2) < 1 ?
      R[j]/sb2 :
      sign(t2)*(delta1[j]-delta2[j])*mp::sqrt(nu/t2t2+1)/(t1/t2-1);
    dabminusRoversb2[j] = delta2[j]*asb2 - Roversb2;
    H[0][j] = dnorm128(R[j]);
    M1[0][j] = asb1 * dnormdsb1[j] * pnorm128(dabminusRoversb1[j]);
    M2[0][j] = asb2 * dnormdsb2[j] * pnorm128(dabminusRoversb2[j]);
    aRminusdelta1[j] = asb1*Roversb1-delta1[j];
    aRminusdelta2[j] = asb2*Roversb2-delta2[j];
    Lfactor1[j] = ab1 * dnorm128(aRminusdelta1[j]);
    Lfactor2[j] = ab2 * dnorm128(aRminusdelta2[j]);
  }
  if(nu >= 3){
    for(j=0; j<J; j++){
      H[1][j] = xdnormx(R[j]);
      M1[1][j] = delta1[j]*ab1*M1[0][j] + ab1 * dnormdsb1[j] *
                  dnorm128(dabminusRoversb1[j]);
      M2[1][j] = delta2[j]*ab2*M2[0][j] + ab2 * dnormdsb2[j] *
                  dnorm128(dabminusRoversb2[j]);
    }
    if(nu >= 4){
      size_t k;
      if(algo == 1){
        mp::float128 A[n]; A[0] = 1; A[1] = 1;
        mp::float128 L[n-2][J];
        for(j=0; j<J; j++){
          L[0][j] = 0.5*H[1][j];
        }
        for(k=2; k<n; k++){
          A[k] = 1.0/(k*A[k-1]);
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
            M1[k][j] = (k-1.0)/k *
              (A[k-2] * delta1[j] * ab1 * M1[k-1][j] + b1*M1[k-2][j]) + 
                Lfactor1[j]*L[k-2][j];
            M2[k][j] = (k-1.0)/k *
              (A[k-2] * delta2[j] * ab2 * M2[k-1][j] + b2*M2[k-2][j]) + 
                Lfactor2[j]*L[k-2][j];
          }
        }
      }else{ // algo 2
        mp::float128 A[n-1]; A[0] = mp::float128(1);
        mp::float128 halfRR[J]; mp::float128 logR[j];
        for(j=0; j<J; j++){
          mp::float128 Rj = mp::float128(R[j]);
          halfRR[j] = 0.5*Rj*Rj;
          logR[j] = mp::log(Rj);
        }
        for(size_t k=0; k<n-2; k++){
          A[k+1] = 1.0/((k+1)*A[k]); // un de trop
          mp::float128 ldf;
          if(k % 2 == 0){
            ldf = 0.5*(k+2)*log_two128 + m::lgamma(mp::float128(0.5*k+2.0));
          }else{
            ldf = m::lgamma(mp::float128(k+3)) - 0.5*(mp::float128(k+1))*log_two128 - 
                    m::lgamma(mp::float128(0.5*(k+3)));
          }
          mp::float128 r = mp::float128(k+1)/mp::float128(k+2);
          for(j=0; j<J; j++){
            mp::float128 K =
              mp::exp(-ldf + (k+1)*logR[j] - halfRR[j] - log_root_two_pi128);
            H[k+2][j] = K*R[j];
            M1[k+2][j] = r *
              (A[k] * delta1[j] * ab1 * M1[k+1][j] + b1*M1[k][j]) +
                K*Lfactor1[j];
            M2[k+2][j] = r *
              (A[k] * delta2[j] * ab2 * M2[k+1][j] + b2*M2[k][j]) +
                K*Lfactor2[j];
          }
        }
      }
    }
  }
  if(nu % 2 == 0){
    for(j=0; j<J; j++){
      mp::float128 sumH=0; mp::float128 sumM=0; mp::float128 out128;
      for(size_t i=0; i<n; i+=2){
        sumH += H[i][j];
        sumM += M1[i][j] - M2[i][j];
      }
      out128 = root_two_pi128*(sumM - sumH * 
                (pnorm128(aRminusdelta1[j]) - pnorm128(aRminusdelta2[j])));
      out[j] = out128.convert_to<double>();
    }
    return out;
  }else{
    NumericVector C = OwenCDF2_C(nu, t1, t2, delta1, delta2);
    for(j=0; j<J; j++){
      mp::float128 sumH=0; mp::float128 sumM=0; mp::float128 out128;
      for(size_t i=1; i<n; i+=2){
        sumH += H[i][j];
        sumM += M1[i][j] - M2[i][j];
      }
      out128 = 2*(sumM - sumH * 
                (pnorm128(aRminusdelta1[j]) - pnorm128(aRminusdelta2[j])));
      out[j] = C[j] + out128.convert_to<double>();
    }
    return out;
  }
}

// --- Owen cumulative function 1 ------------------------------------------- //
NumericVector OwenCDF1_C
    (int nu, double t1, double t2, NumericVector delta1, NumericVector delta2){
  const double a1 = sign(t1)*sqrt(t1*t1/nu);
  const double sb1 = sqrt(nu/(nu+t1*t1));
  const double sb2 = sqrt(nu/(nu+t2*t2));
  size_t J = delta1.size();
  NumericVector C(J);
  for(size_t j=0; j<J; j++){
    double R = sqrt(nu)*(delta1[j] - delta2[j])/(t1-t2);
    double C1 = owent(delta1[j]*sb1, a1);
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
    C[j] = 2*(C1 + C2 + C3) - (delta1[j] >= 0) + (delta2[j] >= 0) +
            pnorm64(-delta2[j]*sb2);
  }
  return C;
}

// [[Rcpp::export]]
NumericVector OwenCDF1
    (size_t nu, double t1, double t2, NumericVector delta1, NumericVector delta2, int algo=1){
  if(nu == 1){
    return OwenCDF1_C(nu, t1, t2, delta1, delta2);
  }
  const size_t J = delta1.size();
  size_t j;
  NumericVector out(J);
  const mp::float128 t1t1(t1*t1);
  const mp::float128 b1 = nu/(nu+t1t1);
  const mp::float128 sb1 = mp::sqrt(b1);
  const mp::float128 ab1 = mp::float128(sqrt(nu)/(nu/t1+t1));
  const mp::float128 asb1 = sign(t1)/mp::sqrt(nu/t1t1+1);
  const mp::float128 t2t2(t2*t2);
  const mp::float128 b2 = nu/(nu+t2t2);
  const mp::float128 sb2 = mp::sqrt(b2);
  const mp::float128 ab2 = mp::float128(sqrt(nu)/(nu/t2+t2));
  const mp::float128 asb2 = sign(t2)/mp::sqrt(nu/t2t2+1);
  mp::float128 R[J];
  mp::float128 dnormdsb1[J]; mp::float128 dnormdsb2[J];
  mp::float128 dabminusRoversb1[J]; mp::float128 dabminusRoversb2[J];
  const size_t n = nu-1;
  mp::float128 M1[n][J]; mp::float128 M2[n][J];
  mp::float128 Lfactor1[J]; mp::float128 Lfactor2[J];
  for(j=0; j<J; j++){
    R[j] = mp::float128(sqrt(nu)*(delta1[j]-delta2[j])/(t1-t2));
    dnormdsb1[j] = dnorm128(delta1[j] * sb1);
    dnormdsb2[j] = dnorm128(delta2[j] * sb2);
    mp::float128 Roversb1 = fabs(t1) < 1 ?
      R[j]/sb1 :
      sign(t1)*(delta1[j]-delta2[j])*mp::sqrt(nu/t1t1+1)/(1-t2/t1);
    dabminusRoversb1[j] = delta1[j]*asb1 - Roversb1;
    mp::float128 Roversb2 = fabs(t2) < 1 ?
      R[j]/sb2 :
      sign(t2)*(delta1[j]-delta2[j])*mp::sqrt(nu/t2t2+1)/(t1/t2-1);
    dabminusRoversb2[j] = delta2[j]*asb2 - Roversb2;
    M1[0][j] = asb1 * dnormdsb1[j] *
      (pnorm128(delta1[j]*asb1) - pnorm128(dabminusRoversb1[j]));
    M2[0][j] = asb2 * dnormdsb2[j] * pnorm128(dabminusRoversb2[j]);
    Lfactor1[j] = ab1 * dnorm128(asb1*Roversb1-delta1[j]);
    Lfactor2[j] = ab2 * dnorm128(asb2*Roversb2-delta2[j]);
  }
  if(nu >= 3){
    for(j=0; j<J; j++){
      M1[1][j] = delta1[j]*ab1*M1[0][j] + ab1 * dnormdsb1[j] *
          (dnorm128(delta1[j]*asb1) - dnorm128(dabminusRoversb1[j]));
      M2[1][j] = delta2[j]*ab2*M2[0][j] + ab2 * dnormdsb2[j] *
                  dnorm128(dabminusRoversb2[j]);
    }
    if(nu >= 4){
      size_t k;
      if(algo == 1){
        mp::float128 A[n]; A[0] = 1; A[1] = 1;
        mp::float128 L[n-2][J];
        for(j=0; j<J; j++){
          L[0][j] = 0.5*xdnormx(R[j]);
        }
        for(k=2; k<n; k++){
          A[k] = 1.0/(k*A[k-1]);
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
            M1[k][j] = (k-1.0)/k *
              (A[k-2] * delta1[j] * ab1 * M1[k-1][j] + b1*M1[k-2][j]) - 
                Lfactor1[j]*L[k-2][j];
            M2[k][j] = (k-1.0)/k *
              (A[k-2] * delta2[j] * ab2 * M2[k-1][j] + b2*M2[k-2][j]) + 
                Lfactor2[j]*L[k-2][j];
          }
        }
      }else{ // algo 2
        mp::float128 A[n-1]; A[0] = mp::float128(1);
        mp::float128 halfRR[J]; mp::float128 logR[j];
        for(j=0; j<J; j++){
          mp::float128 Rj = mp::float128(R[j]);
          halfRR[j] = 0.5*Rj*Rj;
          logR[j] = mp::log(Rj);
        }
        for(k=0; k<n-2; k++){
          A[k+1] = 1.0/((k+1)*A[k]); // un de trop
          mp::float128 ldf;
          if(k % 2 == 0){
            ldf = 0.5*(k+2)*log_two128 + m::lgamma(mp::float128(0.5*k+2.0));
          }else{
            ldf = m::lgamma(mp::float128(k+3)) - 0.5*(mp::float128(k+1))*log_two128 - 
                    m::lgamma(mp::float128(0.5*(k+3)));
          }
          mp::float128 r = mp::float128(k+1)/mp::float128(k+2);
          for(j=0; j<J; j++){
            mp::float128 K =
              mp::exp(-ldf + (k+1)*logR[j] - halfRR[j] - log_root_two_pi128);
            M1[k+2][j] = r *
              (A[k] * delta1[j] * ab1 * M1[k+1][j] + b1*M1[k][j]) -
                K*Lfactor1[j];
            M2[k+2][j] = r *
              (A[k] * delta2[j] * ab2 * M2[k+1][j] + b2*M2[k][j]) +
                K*Lfactor2[j];
          }
        }
      }
    }
  }
  if(nu % 2 == 0){
    for(j=0; j<J; j++){
      mp::float128 sumM=0; mp::float128 out128;
      for(size_t i=0; i<n; i+=2){
        sumM += M2[i][j] + M1[i][j];
      }
      out128 = root_two_pi128*sumM + pnorm128(-delta1[j]);
      out[j] = out128.convert_to<double>();
    }
    return out;
  }else{
    NumericVector C = OwenCDF1_C(nu, t1, t2, delta1, delta2);
    for(j=0; j<J; j++){
      mp::float128 sumM=0;
      for(size_t i=1; i<n; i+=2){
        sumM += M2[i][j] + M1[i][j];
      }
      out[j] = C[j] + (2*sumM).convert_to<double>();
    }
    return out;
  }
}

// --- Owen cumulative function 3 ------------------------------------------- //
NumericVector OwenCDF3_C
    (int nu, double t1, double t2, NumericVector delta1, NumericVector delta2){
  const double sb1 = sqrt(nu/(nu+t1*t1));
  const double a2 = sign(t2)*sqrt(t2*t2/nu);
  const double sb2 = sqrt(nu/(nu+t2*t2));
  size_t J = delta1.size();
  NumericVector C(J);
  for(size_t j=0; j<J; j++){
    double R = sqrt(nu)*(delta1[j] - delta2[j])/(t1-t2);
    double C1 = -owent(delta2[j]*sb2, a2);
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
    C[j] = 2*(C1 + C2 + C3) - (delta1[j] >= 0) + (delta2[j] >= 0) -
      pnorm64(-delta1[j]*sb1) + 1;
  }
  return C;
}

// [[Rcpp::export]]
NumericVector OwenCDF3
    (size_t nu, double t1, double t2, NumericVector delta1, NumericVector delta2, int algo=1){
  if(nu == 1){
    return OwenCDF3_C(nu, t1, t2, delta1, delta2);
  }
  const size_t J = delta1.size();
  size_t j;
  NumericVector out(J);
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
  mp::float128 dnormdsb1[J]; mp::float128 dnormdsb2[J];
  mp::float128 aRminusdelta1[J]; mp::float128 aRminusdelta2[J];
  mp::float128 dabminusRoversb1[J]; mp::float128 dabminusRoversb2[J];
  const size_t n = nu-1;
  mp::float128 H[n][J]; mp::float128 M1[n][J]; mp::float128 M2[n][J];
  mp::float128 Lfactor1[J]; mp::float128 Lfactor2[J];
  for(j=0; j<J; j++){
    R[j] = mp::float128(sqrt(nu)*(delta1[j]-delta2[j])/(t1-t2));
    dnormdsb1[j] = dnorm128(delta1[j] * sb1);
    dnormdsb2[j] = dnorm128(delta2[j] * sb2);
    mp::float128 Roversb1 = fabs(t1) < 1 ?
      R[j]/sb1 :
      sign(t1)*(delta1[j]-delta2[j])*mp::sqrt(nu/t1t1+1)/(1-t2/t1);
    dabminusRoversb1[j] = delta1[j]*asb1 - Roversb1;
    mp::float128 Roversb2 = fabs(t2) < 1 ?
      R[j]/sb2 :
      sign(t2)*(delta1[j]-delta2[j])*mp::sqrt(nu/t2t2+1)/(t1/t2-1);
    dabminusRoversb2[j] = delta2[j]*asb2 - Roversb2;
    H[0][j] = dnorm128(R[j]);
    M1[0][j] = asb1 * dnormdsb1[j] * pnorm128(dabminusRoversb1[j]);
    M2[0][j] = asb2 * dnormdsb2[j] *
      (pnorm128(delta2[j]*asb2) - pnorm128(dabminusRoversb2[j]));
    aRminusdelta1[j] = asb1*Roversb1-delta1[j];
    aRminusdelta2[j] = asb2*Roversb2-delta2[j];
    Lfactor1[j] = ab1 * dnorm128(aRminusdelta1[j]);
    Lfactor2[j] = ab2 * dnorm128(aRminusdelta2[j]);
  }
  if(nu >= 3){
    for(j=0; j<J; j++){
      H[1][j] = xdnormx(R[j]);
      M1[1][j] = delta1[j]*ab1*M1[0][j] + ab1 * dnormdsb1[j] *
                  dnorm128(dabminusRoversb1[j]);
      M2[1][j] = delta2[j]*ab2*M2[0][j] + ab2 * dnormdsb2[j] *
                  (dnorm128(delta2[j]*asb2) - dnorm128(dabminusRoversb2[j]));
    }
    if(nu >= 4){
      size_t k;
      if(algo == 1){
        mp::float128 A[n]; A[0] = 1; A[1] = 1;
        mp::float128 L[n-2][J];
        for(j=0; j<J; j++){
          L[0][j] = 0.5*H[1][j];
        }
        for(k=2; k<n; k++){
          A[k] = 1.0/(k*A[k-1]);
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
            M1[k][j] = (k-1.0)/k *
              (A[k-2] * delta1[j] * ab1 * M1[k-1][j] + b1*M1[k-2][j]) + 
                Lfactor1[j]*L[k-2][j];
            M2[k][j] = (k-1.0)/k *
              (A[k-2] * delta2[j] * ab2 * M2[k-1][j] + b2*M2[k-2][j]) - 
                Lfactor2[j]*L[k-2][j];
          }
        }
      }else{
        mp::float128 A[n-1]; A[0] = mp::float128(1);
        mp::float128 halfRR[J]; mp::float128 logR[j];
        for(j=0; j<J; j++){
          mp::float128 Rj = mp::float128(R[j]);
          halfRR[j] = 0.5*Rj*Rj;
          logR[j] = mp::log(Rj);
        }
        for(k=0; k<n-2; k++){
          A[k+1] = 1.0/((k+1)*A[k]); // un de trop
          mp::float128 ldf;
          if(k % 2 == 0){
            ldf = 0.5*(k+2)*log_two128 + m::lgamma(mp::float128(0.5*k+2.0));
          }else{
            ldf = m::lgamma(mp::float128(k+3)) - 0.5*(mp::float128(k+1))*log_two128 - 
                    m::lgamma(mp::float128(0.5*(k+3)));
          }
          mp::float128 r = mp::float128(k+1)/mp::float128(k+2);
          for(j=0; j<J; j++){
            mp::float128 K =
              mp::exp(-ldf + (k+1)*logR[j] - halfRR[j] - log_root_two_pi128);
            M1[k+2][j] = r *
              (A[k] * delta1[j] * ab1 * M1[k+1][j] + b1*M1[k][j]) +
                K*Lfactor1[j];
            M2[k+2][j] = r *
              (A[k] * delta2[j] * ab2 * M2[k+1][j] + b2*M2[k][j]) -
                K*Lfactor2[j];
          }
        }
      }
    }
  }
  if(nu % 2 == 0){
    for(j=0; j<J; j++){
      mp::float128 sumH=0; mp::float128 sumM=0; mp::float128 out128;
      for(size_t i=0; i<n; i+=2){
        sumH += H[i][j];
        sumM += -M1[i][j] - M2[i][j];
      }
      out128 = root_two_pi128*(sumM + sumH *
                (pnorm128(aRminusdelta1[j]) - pnorm128(aRminusdelta2[j]))) +
                1 - pnorm128(-delta2[j]);
      out[j] = out128.convert_to<double>();
    }
    return out;
  }else{
    NumericVector C = OwenCDF3_C(nu, t1, t2, delta1, delta2);
    for(j=0; j<J; j++){
      mp::float128 sumH=0; mp::float128 sumM=0; mp::float128 out128;
      for(size_t i=1; i<n; i+=2){
        sumH += H[i][j];
        sumM += -M1[i][j] - M2[i][j];
      }
      out128 = 2*(sumM + sumH * 
                (pnorm128(aRminusdelta1[j]) - pnorm128(aRminusdelta2[j])));
      out[j] = C[j] + out128.convert_to<double>();
    }
    return out;
  }
}
