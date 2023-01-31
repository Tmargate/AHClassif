#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix dissim(NumericMatrix X, String method = "average", bool dissim_matrix = false){
  if((method != "single") && (method != "complete") && (method != "average") && (method != "weighted") && (method != "centroid") && (method != "median") && (method != "ward")){
    return false;
  }
  int n0 = X.nrow();
  NumericMatrix D(n0,n0);
  if(dissim_matrix == false){
    for(int j=0; j< (n0-1); j++){
      for(int i=j; i < n0; i++) {
        NumericVector v = X(i,_) - X(j,_);
        double sum =0;
        for(int k=0; k<n0;k++){
          sum = sum + v[i]*v[i];
        }
        D(i,j) = sum;
        D(j,i) = sum;
      }
    }
  }
  else {
    D = X;
    for(int i=0; i<n0;i++){
      for(int j =0; j< n0; j++){
        if(D(i,j) != D(j,i)){
          return false;
        }
      }
    }
  }
  /* manque le inf */
  double alpha_i = 0.5;
  double alpha_j = 0.5;
  double beta = 0;
  double gamma = 0;
  if(method == "single"){
    gamma =0.5;
  }
  if(method == "complete"){
    gamma = 0.5;
  }
  if(method == "median"){
    beta = -0.25;
  }
  NumericMatrix P(n0,n0);
  for(int i=0; i< n0; i++){
    P(1,i) = i+1;
  }

  int n = n0;
  for(int k =1; k<n; k++){
    int mini = D(0,0);
    int argmin_i = 0;
    int argmin_j = 0;
    for(int p =0; p<n;p++){
      for(int q =0; q<n;q++){
        if(D(q,p) < mini){
          argmin_i = q+1;
          argmin_j = p+1;
          mini=D(q,p);
        }
      }
    }
    int j = argmin_j;
    int i = argmin_i;
    P(k,_) = P(k+1,_);
    for(int p =0; p<n; p++){
      if(P(k,p)==i){
        P(k,p)=j;
      }
      if(P(k,p) > i){
        P(k,p) = P(k,p) -1;
      }
    }
    if((method =="average") || (method == "centroid") || (method == "ward")){
      int sum_i = 0;
      int sum_j = 0;
      for(int p =0; p<n;p++){
        if(P(k-1,p) == i){
          sum_i = sum_i +1;
        }
        if(P(k-1,p)==j){
          sum_j = sum_j +1;
        }
      }
      int card_i = sum_i;
      int card_j = sum_j;
      if((method == "average") || (method == "centroid")){
        alpha_i = card_i/(card_i + card_j);
        alpha_j = card_j/(card_i + card_j);
        if(method=="centroid"){
          beta = -alpha_i*alpha_j;
        }
      }
      else{
        NumericVector alpha_i(n);
        NumericVector alpha_j(n);
        NumericVector beta(n);
        for(int p =0; p<n;p++){
          alpha_i[p] = card_i;
          alpha_j[p] = card_j;
          beta[p] = 0;
        }
        for(int m=0;m<n;m++){
          int card_m = 0;
          for(int p =0; p<n; p++){
            if(P(k-1,p)==m){
              card_m = card_m +1;
            }
          }
          int S = card_i + card_j + card_m;
          alpha_i[m] = (alpha_i[m] + card_m)/S;
          alpha_j[m] = (alpha_j[m] + card_m)/S;
          beta[m] = -card_m/S;
        }
      }
    }
    for(int p =0;p<n;p++){
      D(j,p) = alpha_i*D(i,p)+alpha_j*D(j,p)+beta*D(i,j)+gamma*abs(D(i,p) - D(j,p));
      /* mettre des infs */
      D(p,j) = D(j,p);
      /* D(-1,-1) ? */
      n = n-1;
    }
  }
  return P;
}
