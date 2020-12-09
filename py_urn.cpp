// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat PY_urn(int n, double alpha, double sigma){
  
  // Frequency start at 1
  arma::vec freq(1); freq(0) = 1;
  
  // Distinct values among x
  int k; k = freq.n_elem;
  
  // Estimate the alpha parameter
  for (int i = 1; i < n; i++){
    
    // Sample X_i | X_1,...,X_{i-1}    
    double new_x; 
    new_x = R::rbinom(1, (alpha + k*sigma)/(alpha + i)); // It's i and not i-1 because of the C++ indexing
    
    // Suppose that X_i is a new value
    if(new_x==1) {
      
      freq  = join_vert(freq,arma::ones(1)); // Add an element to the vector
      k = freq.n_elem; // Increase the number of distinct values
    } else {
      
      // Suppose instead X_i = X^*_j, for some j = 1,...,k
      Rcpp::IntegerVector clusters = Rcpp::Range(1,k); // Range of possibile clusters
      
      // Which is the X^*_j?
      arma::vec probs; probs = (freq - sigma)/sum(freq - sigma);
      int idx; idx = RcppArmadillo::sample(clusters, 1, TRUE, probs)[0]-1; 
      // The frequency of that class is increased
      freq[idx] = freq[idx] + 1;
    }
  }
  //out = join_horiz(distinct_species,families);
  return(freq);
}

