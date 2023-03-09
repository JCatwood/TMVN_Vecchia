#ifndef MVPHI_H
#define MVPHI_H
#include <Rcpp.h>
#include <R.h>
#include <Rinternals.h>
#include <algorithm>

#ifdef __cplusplus
extern "C" {
#endif

void F77_NAME(mvphnv)(double *p , double *v);

void F77_NAME(mvphi)(double *z , double *v);

#ifdef __cplusplus
}
#endif

inline void lc_vdCdfNorm(int N , double *z , double *v) {
	for(int i = 0 ; i < N ; i++) F77_CALL(mvphi)(z+i , v+i);
}

inline void lc_vdCdfNormInv(int N , double *p , double *v) {
	for(int i = 0 ; i < N ; i++) F77_CALL(mvphnv)(p+i , v+i);
}

inline void lc_vdCdfNorm(int N , double *z , double *v , int stride) {
	for(int i = 0 ; i < N ; i++) {
		F77_CALL(mvphi)(z , v);
		z += stride;
		v += stride;
	}
}

inline void lc_vdCdfNormInv(int N , double *p , double *v , int stride) {
	for(int i = 0 ; i < N ; i++) {
		F77_CALL(mvphnv)(p , v);
		p += stride;
		v += stride;
	}
}

inline void lc_vdCdfNormInv(int N , double *p , int stride_p, double *v , int stride_v) {
	for(int i = 0 ; i < N ; i++) {
		F77_CALL(mvphnv)(p , v);
		p += stride_p;
		v += stride_v;
	}
}
#endif
