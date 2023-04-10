#include <Rcpp.h>
#include <algorithm>


using namespace Rcpp;


/*
    Find ordered nearest neighbors based on a correlation Matrix. 
    Assuming the absolute value of the correlation is monotonically 
    decreasing with distance.
    Returns an n X (m + 1) matrix similar to `GpGp::find_ordered_nn`.
*/
// [[Rcpp::export]]
IntegerMatrix find_nn_corr_internal(const NumericMatrix &corrMat, int m){
	int n = corrMat.rows();
	IntegerMatrix NN(n, m + 1);
	NN.fill(NA_INTEGER);
	NN(0, 0) = 0;
	for(int i = 1; i < n; i++){
		int *order = new int[i + 1];
		std::iota(order, order + i + 1, 0);
		std::(order, order + i + 1, [&corrMat, &i](int &j, int &k){
			return corrMat(j, i) > corrMat(k, i); });
		for(int j = 0; j < min(m + 1, i + 1); j++)
			NN(i, j) = order[j];
		delete[] order;
	}
	return NN;
}
