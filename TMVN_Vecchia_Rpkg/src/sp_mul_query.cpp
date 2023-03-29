#include <Rcpp.h>


using namespace Rcpp;
using namespace std;


/*
    Query entries of a sparse matrix multiplied by its conjuate.
    Given sparse matrix M, for each query index pair (query_row, query_col),
      if M is CSR:
        compute M[query_row, ] * M[query_col, ], equivalent to
          (M \cdot M^{\top})[query_row, query_col]
      if M is CSC:
        compute M[, query_row] * M[, query_col], equivalent to
          (M^{\top} \cdot M)[query_row, query_col]
    Implementations of both cases are the same
*/
// [[Rcpp::export]]
NumericVector sp_mat_mul_query(
	const IntegerVector &queryRow, const IntegerVector &queryCol,
    const IntegerVector &idx, const IntegerVector &cidx,
    const NumericVector &val)
{
	int n_query = queryRow.size();
	NumericVector query_val(n_query);
	// parallelization can be applied here
	for(int k = 0; k < n_query; k++){
		int i = queryRow[k] - 1;
		int j = queryCol[k] - 1;
		int i_bgn = cidx[i];
		int i_end = cidx[i + 1];
		int j_bgn = cidx[j];
		int j_end = cidx[j + 1];
		double sum = 0;
		while((i_bgn < i_end) && (j_bgn < j_end)){
			if(idx[i_bgn] == idx[j_bgn]){
				sum +=  val[i_bgn] * val[j_bgn];
				i_bgn++;
				j_bgn++;
			}else if(idx[i_bgn] < idx[j_bgn]){
				i_bgn++;
			}else{
				j_bgn++;
			}
		}
		query_val[k] = sum;
	}
	return query_val;
}
