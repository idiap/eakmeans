
#ifdef WITHBLAS

#include "blastemplates.h"
#include <cblas.h>


namespace wblas{
template <> 
float dot(OPENBLAS_CONST blasint n, OPENBLAS_CONST float *x, OPENBLAS_CONST blasint incx, OPENBLAS_CONST float *y, OPENBLAS_CONST blasint incy){
	return cblas_sdot(n, x, incx, y, incy);
}

template <> 
double dot(OPENBLAS_CONST blasint n, OPENBLAS_CONST double *x, OPENBLAS_CONST blasint incx, OPENBLAS_CONST double *y, OPENBLAS_CONST blasint incy){
	return cblas_ddot(n, x, incx, y, incy);
}

template <> 
float asum(OPENBLAS_CONST blasint n, OPENBLAS_CONST float *x, OPENBLAS_CONST blasint incx){
	return cblas_sasum(n, x, incx);
}

template <> 
double asum(OPENBLAS_CONST blasint n, OPENBLAS_CONST double *x, OPENBLAS_CONST blasint incx){
	return cblas_dasum(n, x, incx);
}

template <> 
float nrm2(OPENBLAS_CONST blasint N, OPENBLAS_CONST float *X, OPENBLAS_CONST blasint incX){
	return cblas_snrm2(N, X, incX);
}

template <> 
double nrm2(OPENBLAS_CONST blasint N, OPENBLAS_CONST double *X, OPENBLAS_CONST blasint incX){
	return cblas_dnrm2(N, X, incX);
}

template <> 
CBLAS_INDEX i_amax(OPENBLAS_CONST blasint n, OPENBLAS_CONST float *x, OPENBLAS_CONST blasint incx){
	return cblas_isamax(n, x, incx);
}

template <> 
CBLAS_INDEX i_amax(OPENBLAS_CONST blasint n, OPENBLAS_CONST double *x, OPENBLAS_CONST blasint incx){
	return cblas_idamax(n, x, incx);
}

template <> 
void axpy(OPENBLAS_CONST blasint n, OPENBLAS_CONST float alpha, OPENBLAS_CONST float *x, OPENBLAS_CONST blasint incx, float *y, OPENBLAS_CONST blasint incy){
	cblas_saxpy(n, alpha, x, incx, y, incy);
}

template <> 
void axpy(OPENBLAS_CONST blasint n, OPENBLAS_CONST double alpha, OPENBLAS_CONST double *x, OPENBLAS_CONST blasint incx, double *y, OPENBLAS_CONST blasint incy){
	cblas_daxpy(n, alpha, x, incx, y, incy);
}

template <> 
void copy(OPENBLAS_CONST blasint n, OPENBLAS_CONST float *x, OPENBLAS_CONST blasint incx, float *y, OPENBLAS_CONST blasint incy){
	cblas_scopy(n, x, incx, y, incy);
}

template <> 
void copy(OPENBLAS_CONST blasint n, OPENBLAS_CONST double *x, OPENBLAS_CONST blasint incx, double *y, OPENBLAS_CONST blasint incy){
	cblas_dcopy(n, x, incx, y, incy);
}

template <> 
void swap(OPENBLAS_CONST blasint n, float *x, OPENBLAS_CONST blasint incx, float *y, OPENBLAS_CONST blasint incy){
	cblas_sswap(n, x, incx, y, incy);
}

template <> 
void swap(OPENBLAS_CONST blasint n, double *x, OPENBLAS_CONST blasint incx, double *y, OPENBLAS_CONST blasint incy){
	cblas_dswap(n, x, incx, y, incy);
}

template <> 
void rot(OPENBLAS_CONST blasint N, float *X, OPENBLAS_CONST blasint incX, float *Y, OPENBLAS_CONST blasint incY, OPENBLAS_CONST float c, OPENBLAS_CONST float s){
	cblas_srot(N, X, incX, Y, incY, c, s);
}

template <> 
void rot(OPENBLAS_CONST blasint N, double *X, OPENBLAS_CONST blasint incX, double *Y, OPENBLAS_CONST blasint incY, OPENBLAS_CONST double c, OPENBLAS_CONST double s){
	cblas_drot(N, X, incX, Y, incY, c, s);
}

template <> 
void rotg(float *a, float *b, float *c, float *s){
	cblas_srotg(a, b, c, s);
}

template <> 
void rotg(double *a, double *b, double *c, double *s){
	cblas_drotg(a, b, c, s);
}

template <> 
void rotm(OPENBLAS_CONST blasint N, float *X, OPENBLAS_CONST blasint incX, float *Y, OPENBLAS_CONST blasint incY, OPENBLAS_CONST float *P){
	cblas_srotm(N, X, incX, Y, incY, P);
}

template <> 
void rotm(OPENBLAS_CONST blasint N, double *X, OPENBLAS_CONST blasint incX, double *Y, OPENBLAS_CONST blasint incY, OPENBLAS_CONST double *P){
	cblas_drotm(N, X, incX, Y, incY, P);
}

template <> 
void rotmg(float *d1, float *d2, float *b1, OPENBLAS_CONST float b2, float *P){
	cblas_srotmg(d1, d2, b1, b2, P);
}

template <> 
void rotmg(double *d1, double *d2, double *b1, OPENBLAS_CONST double b2, double *P){
	cblas_drotmg(d1, d2, b1, b2, P);
}

template <> 
void scal(OPENBLAS_CONST blasint N, OPENBLAS_CONST float alpha, float *X, OPENBLAS_CONST blasint incX){
	cblas_sscal(N, alpha, X, incX);
}

template <> 
void scal(OPENBLAS_CONST blasint N, OPENBLAS_CONST double alpha, double *X, OPENBLAS_CONST blasint incX){
	cblas_dscal(N, alpha, X, incX);
}

template <> 
void gemv(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_TRANSPOSE trans, OPENBLAS_CONST blasint m, OPENBLAS_CONST blasint n,OPENBLAS_CONST float alpha, OPENBLAS_CONST float *a, OPENBLAS_CONST blasint lda, OPENBLAS_CONST float *x, OPENBLAS_CONST blasint incx, OPENBLAS_CONST float beta, float *y, OPENBLAS_CONST blasint incy){
	cblas_sgemv(order, trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
}

template <> 
void gemv(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_TRANSPOSE trans, OPENBLAS_CONST blasint m, OPENBLAS_CONST blasint n,OPENBLAS_CONST double alpha, OPENBLAS_CONST double *a, OPENBLAS_CONST blasint lda, OPENBLAS_CONST double *x, OPENBLAS_CONST blasint incx, OPENBLAS_CONST double beta, double *y, OPENBLAS_CONST blasint incy){
	cblas_dgemv(order, trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
}

template <> 
void ger(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST blasint M, OPENBLAS_CONST blasint N, OPENBLAS_CONST float alpha, OPENBLAS_CONST float *X, OPENBLAS_CONST blasint incX, OPENBLAS_CONST float *Y, OPENBLAS_CONST blasint incY, float *A, OPENBLAS_CONST blasint lda){
	cblas_sger(order, M, N, alpha, X, incX, Y, incY, A, lda);
}

template <> 
void ger(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST blasint M, OPENBLAS_CONST blasint N, OPENBLAS_CONST double alpha, OPENBLAS_CONST double *X, OPENBLAS_CONST blasint incX, OPENBLAS_CONST double *Y, OPENBLAS_CONST blasint incY, double *A, OPENBLAS_CONST blasint lda){
	cblas_dger(order, M, N, alpha, X, incX, Y, incY, A, lda);
}

template <> 
void trsv(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA, OPENBLAS_CONST enum CBLAS_DIAG Diag, OPENBLAS_CONST blasint N, OPENBLAS_CONST float *A, OPENBLAS_CONST blasint lda, float *X, OPENBLAS_CONST blasint incX){
	cblas_strsv(order, Uplo, TransA, Diag, N, A, lda, X, incX);
}

template <> 
void trsv(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA, OPENBLAS_CONST enum CBLAS_DIAG Diag, OPENBLAS_CONST blasint N, OPENBLAS_CONST double *A, OPENBLAS_CONST blasint lda, double *X, OPENBLAS_CONST blasint incX){
	cblas_dtrsv(order, Uplo, TransA, Diag, N, A, lda, X, incX);
}

template <> 
void trmv(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA, OPENBLAS_CONST enum CBLAS_DIAG Diag, OPENBLAS_CONST blasint N, OPENBLAS_CONST float *A, OPENBLAS_CONST blasint lda, float *X, OPENBLAS_CONST blasint incX){
	cblas_strmv(order, Uplo, TransA, Diag, N, A, lda, X, incX);
}

template <> 
void trmv(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA, OPENBLAS_CONST enum CBLAS_DIAG Diag, OPENBLAS_CONST blasint N, OPENBLAS_CONST double *A, OPENBLAS_CONST blasint lda, double *X, OPENBLAS_CONST blasint incX){
	cblas_dtrmv(order, Uplo, TransA, Diag, N, A, lda, X, incX);
}

template <> 
void syr(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST blasint N, OPENBLAS_CONST float alpha, OPENBLAS_CONST float *X, OPENBLAS_CONST blasint incX, float *A, OPENBLAS_CONST blasint lda){
	cblas_ssyr(order, Uplo, N, alpha, X, incX, A, lda);
}

template <> 
void syr(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST blasint N, OPENBLAS_CONST double alpha, OPENBLAS_CONST double *X, OPENBLAS_CONST blasint incX, double *A, OPENBLAS_CONST blasint lda){
	cblas_dsyr(order, Uplo, N, alpha, X, incX, A, lda);
}

template <> 
void syr2(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_UPLO Uplo,OPENBLAS_CONST blasint N, OPENBLAS_CONST float alpha, OPENBLAS_CONST float *X,OPENBLAS_CONST blasint incX, OPENBLAS_CONST float *Y, OPENBLAS_CONST blasint incY, float *A, OPENBLAS_CONST blasint lda){
	cblas_ssyr2(order, Uplo, N, alpha, X, incX, Y, incY, A, lda);
}

template <> 
void syr2(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST blasint N, OPENBLAS_CONST double alpha, OPENBLAS_CONST double *X,OPENBLAS_CONST blasint incX, OPENBLAS_CONST double *Y, OPENBLAS_CONST blasint incY, double *A, OPENBLAS_CONST blasint lda){
	cblas_dsyr2(order, Uplo, N, alpha, X, incX, Y, incY, A, lda);
}

template <> 
void gbmv(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA, OPENBLAS_CONST blasint M, OPENBLAS_CONST blasint N,OPENBLAS_CONST blasint KL, OPENBLAS_CONST blasint KU, OPENBLAS_CONST float alpha, OPENBLAS_CONST float *A, OPENBLAS_CONST blasint lda, OPENBLAS_CONST float *X, OPENBLAS_CONST blasint incX, OPENBLAS_CONST float beta, float *Y, OPENBLAS_CONST blasint incY){
	cblas_sgbmv(order, TransA, M, N, KL, KU, alpha, A, lda, X, incX, beta, Y, incY);
}

template <> 
void gbmv(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA, OPENBLAS_CONST blasint M, OPENBLAS_CONST blasint N,OPENBLAS_CONST blasint KL, OPENBLAS_CONST blasint KU, OPENBLAS_CONST double alpha, OPENBLAS_CONST double *A, OPENBLAS_CONST blasint lda, OPENBLAS_CONST double *X, OPENBLAS_CONST blasint incX, OPENBLAS_CONST double beta, double *Y, OPENBLAS_CONST blasint incY){
	cblas_dgbmv(order, TransA, M, N, KL, KU, alpha, A, lda, X, incX, beta, Y, incY);
}

template <> 
void sbmv(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST blasint N, OPENBLAS_CONST blasint K, OPENBLAS_CONST float alpha, OPENBLAS_CONST float *A,OPENBLAS_CONST blasint lda, OPENBLAS_CONST float *X, OPENBLAS_CONST blasint incX, OPENBLAS_CONST float beta, float *Y, OPENBLAS_CONST blasint incY){
	cblas_ssbmv(order, Uplo, N, K, alpha, A, lda, X, incX, beta, Y, incY);
}

template <> 
void sbmv(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST blasint N, OPENBLAS_CONST blasint K, OPENBLAS_CONST double alpha, OPENBLAS_CONST double *A,OPENBLAS_CONST blasint lda, OPENBLAS_CONST double *X, OPENBLAS_CONST blasint incX, OPENBLAS_CONST double beta, double *Y, OPENBLAS_CONST blasint incY){
	cblas_dsbmv(order, Uplo, N, K, alpha, A, lda, X, incX, beta, Y, incY);
}

template <> 
void tbmv(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA, OPENBLAS_CONST enum CBLAS_DIAG Diag,OPENBLAS_CONST blasint N, OPENBLAS_CONST blasint K, OPENBLAS_CONST float *A, OPENBLAS_CONST blasint lda, float *X, OPENBLAS_CONST blasint incX){
	cblas_stbmv(order, Uplo, TransA, Diag, N, K, A, lda, X, incX);
}

template <> 
void tbmv(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA, OPENBLAS_CONST enum CBLAS_DIAG Diag,OPENBLAS_CONST blasint N, OPENBLAS_CONST blasint K, OPENBLAS_CONST double *A, OPENBLAS_CONST blasint lda, double *X, OPENBLAS_CONST blasint incX){
	cblas_dtbmv(order, Uplo, TransA, Diag, N, K, A, lda, X, incX);
}

template <> 
void tbsv(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA, OPENBLAS_CONST enum CBLAS_DIAG Diag,OPENBLAS_CONST blasint N, OPENBLAS_CONST blasint K, OPENBLAS_CONST float *A, OPENBLAS_CONST blasint lda, float *X, OPENBLAS_CONST blasint incX){
	cblas_stbsv(order, Uplo, TransA, Diag, N, K, A, lda, X, incX);
}

template <> 
void tbsv(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA, OPENBLAS_CONST enum CBLAS_DIAG Diag,OPENBLAS_CONST blasint N, OPENBLAS_CONST blasint K, OPENBLAS_CONST double *A, OPENBLAS_CONST blasint lda, double *X, OPENBLAS_CONST blasint incX){
	cblas_dtbsv(order, Uplo, TransA, Diag, N, K, A, lda, X, incX);
}

template <> 
void tpmv(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA, OPENBLAS_CONST enum CBLAS_DIAG Diag,OPENBLAS_CONST blasint N, OPENBLAS_CONST float *Ap, float *X, OPENBLAS_CONST blasint incX){
	cblas_stpmv(order, Uplo, TransA, Diag, N, Ap, X, incX);
}

template <> 
void tpmv(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA, OPENBLAS_CONST enum CBLAS_DIAG Diag,OPENBLAS_CONST blasint N, OPENBLAS_CONST double *Ap, double *X, OPENBLAS_CONST blasint incX){
	cblas_dtpmv(order, Uplo, TransA, Diag, N, Ap, X, incX);
}

template <> 
void tpsv(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA, OPENBLAS_CONST enum CBLAS_DIAG Diag,OPENBLAS_CONST blasint N, OPENBLAS_CONST float *Ap, float *X, OPENBLAS_CONST blasint incX){
	cblas_stpsv(order, Uplo, TransA, Diag, N, Ap, X, incX);
}

template <> 
void tpsv(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA, OPENBLAS_CONST enum CBLAS_DIAG Diag,OPENBLAS_CONST blasint N, OPENBLAS_CONST double *Ap, double *X, OPENBLAS_CONST blasint incX){
	cblas_dtpsv(order, Uplo, TransA, Diag, N, Ap, X, incX);
}

template <> 
void symv(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST blasint N, OPENBLAS_CONST float alpha, OPENBLAS_CONST float *A,OPENBLAS_CONST blasint lda, OPENBLAS_CONST float *X, OPENBLAS_CONST blasint incX, OPENBLAS_CONST float beta, float *Y, OPENBLAS_CONST blasint incY){
	cblas_ssymv(order, Uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
}

template <> 
void symv(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST blasint N, OPENBLAS_CONST double alpha, OPENBLAS_CONST double *A,OPENBLAS_CONST blasint lda, OPENBLAS_CONST double *X, OPENBLAS_CONST blasint incX, OPENBLAS_CONST double beta, double *Y, OPENBLAS_CONST blasint incY){
	cblas_dsymv(order, Uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
}

template <> 
void spmv(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST blasint N, OPENBLAS_CONST float alpha, OPENBLAS_CONST float *Ap,OPENBLAS_CONST float *X, OPENBLAS_CONST blasint incX, OPENBLAS_CONST float beta, float *Y, OPENBLAS_CONST blasint incY){
	cblas_sspmv(order, Uplo, N, alpha, Ap, X, incX, beta, Y, incY);
}

template <> 
void spmv(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST blasint N, OPENBLAS_CONST double alpha, OPENBLAS_CONST double *Ap,OPENBLAS_CONST double *X, OPENBLAS_CONST blasint incX, OPENBLAS_CONST double beta, double *Y, OPENBLAS_CONST blasint incY){
	cblas_dspmv(order, Uplo, N, alpha, Ap, X, incX, beta, Y, incY);
}

template <> 
void spr(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST blasint N, OPENBLAS_CONST float alpha, OPENBLAS_CONST float *X, OPENBLAS_CONST blasint incX, float *Ap){
	cblas_sspr(order, Uplo, N, alpha, X, incX, Ap);
}

template <> 
void spr(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST blasint N, OPENBLAS_CONST double alpha, OPENBLAS_CONST double *X, OPENBLAS_CONST blasint incX, double *Ap){
	cblas_dspr(order, Uplo, N, alpha, X, incX, Ap);
}

template <> 
void spr2(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST blasint N, OPENBLAS_CONST float alpha, OPENBLAS_CONST float *X, OPENBLAS_CONST blasint incX, OPENBLAS_CONST float *Y, OPENBLAS_CONST blasint incY, float *A){
	cblas_sspr2(order, Uplo, N, alpha, X, incX, Y, incY, A);
}

template <> 
void spr2(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST blasint N, OPENBLAS_CONST double alpha, OPENBLAS_CONST double *X, OPENBLAS_CONST blasint incX, OPENBLAS_CONST double *Y, OPENBLAS_CONST blasint incY, double *A){
	cblas_dspr2(order, Uplo, N, alpha, X, incX, Y, incY, A);
}

template <> 
void gemm(OPENBLAS_CONST enum CBLAS_ORDER Order, OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA, OPENBLAS_CONST enum CBLAS_TRANSPOSE TransB, OPENBLAS_CONST blasint M, OPENBLAS_CONST blasint N, OPENBLAS_CONST blasint K,OPENBLAS_CONST float alpha, OPENBLAS_CONST float *A, OPENBLAS_CONST blasint lda, OPENBLAS_CONST float *B, OPENBLAS_CONST blasint ldb, OPENBLAS_CONST float beta, float *C, OPENBLAS_CONST blasint ldc){
	cblas_sgemm(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}

template <> 
void gemm(OPENBLAS_CONST enum CBLAS_ORDER Order, OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA, OPENBLAS_CONST enum CBLAS_TRANSPOSE TransB, OPENBLAS_CONST blasint M, OPENBLAS_CONST blasint N, OPENBLAS_CONST blasint K,OPENBLAS_CONST double alpha, OPENBLAS_CONST double *A, OPENBLAS_CONST blasint lda, OPENBLAS_CONST double *B, OPENBLAS_CONST blasint ldb, OPENBLAS_CONST double beta, double *C, OPENBLAS_CONST blasint ldc){
	cblas_dgemm(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}

template <> 
void symm(OPENBLAS_CONST enum CBLAS_ORDER Order, OPENBLAS_CONST enum CBLAS_SIDE Side, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST blasint M, OPENBLAS_CONST blasint N,OPENBLAS_CONST float alpha, OPENBLAS_CONST float *A, OPENBLAS_CONST blasint lda, OPENBLAS_CONST float *B, OPENBLAS_CONST blasint ldb, OPENBLAS_CONST float beta, float *C, OPENBLAS_CONST blasint ldc){
	cblas_ssymm(Order, Side, Uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
}

template <> 
void symm(OPENBLAS_CONST enum CBLAS_ORDER Order, OPENBLAS_CONST enum CBLAS_SIDE Side, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST blasint M, OPENBLAS_CONST blasint N,OPENBLAS_CONST double alpha, OPENBLAS_CONST double *A, OPENBLAS_CONST blasint lda, OPENBLAS_CONST double *B, OPENBLAS_CONST blasint ldb, OPENBLAS_CONST double beta, double *C, OPENBLAS_CONST blasint ldc){
	cblas_dsymm(Order, Side, Uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
}

template <> 
void syrk(OPENBLAS_CONST enum CBLAS_ORDER Order, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST enum CBLAS_TRANSPOSE Trans,OPENBLAS_CONST blasint N, OPENBLAS_CONST blasint K, OPENBLAS_CONST float alpha, OPENBLAS_CONST float *A, OPENBLAS_CONST blasint lda, OPENBLAS_CONST float beta, float *C, OPENBLAS_CONST blasint ldc){
	cblas_ssyrk(Order, Uplo, Trans, N, K, alpha, A, lda, beta, C, ldc);
}

template <> 
void syrk(OPENBLAS_CONST enum CBLAS_ORDER Order, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST enum CBLAS_TRANSPOSE Trans,OPENBLAS_CONST blasint N, OPENBLAS_CONST blasint K, OPENBLAS_CONST double alpha, OPENBLAS_CONST double *A, OPENBLAS_CONST blasint lda, OPENBLAS_CONST double beta, double *C, OPENBLAS_CONST blasint ldc){
	cblas_dsyrk(Order, Uplo, Trans, N, K, alpha, A, lda, beta, C, ldc);
}

template <> 
void syr2k(OPENBLAS_CONST enum CBLAS_ORDER Order, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST enum CBLAS_TRANSPOSE Trans,OPENBLAS_CONST blasint N, OPENBLAS_CONST blasint K, OPENBLAS_CONST float alpha, OPENBLAS_CONST float *A, OPENBLAS_CONST blasint lda, OPENBLAS_CONST float *B, OPENBLAS_CONST blasint ldb, OPENBLAS_CONST float beta, float *C, OPENBLAS_CONST blasint ldc){
	cblas_ssyr2k(Order, Uplo, Trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}

template <> 
void syr2k(OPENBLAS_CONST enum CBLAS_ORDER Order, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST enum CBLAS_TRANSPOSE Trans,OPENBLAS_CONST blasint N, OPENBLAS_CONST blasint K, OPENBLAS_CONST double alpha, OPENBLAS_CONST double *A, OPENBLAS_CONST blasint lda, OPENBLAS_CONST double *B, OPENBLAS_CONST blasint ldb, OPENBLAS_CONST double beta, double *C, OPENBLAS_CONST blasint ldc){
	cblas_dsyr2k(Order, Uplo, Trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}

template <> 
void trmm(OPENBLAS_CONST enum CBLAS_ORDER Order, OPENBLAS_CONST enum CBLAS_SIDE Side, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,OPENBLAS_CONST enum CBLAS_DIAG Diag, OPENBLAS_CONST blasint M, OPENBLAS_CONST blasint N, OPENBLAS_CONST float alpha, OPENBLAS_CONST float *A, OPENBLAS_CONST blasint lda, float *B, OPENBLAS_CONST blasint ldb){
	cblas_strmm(Order, Side, Uplo, TransA, Diag, M, N, alpha, A, lda, B, ldb);
}

template <> 
void trmm(OPENBLAS_CONST enum CBLAS_ORDER Order, OPENBLAS_CONST enum CBLAS_SIDE Side, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,OPENBLAS_CONST enum CBLAS_DIAG Diag, OPENBLAS_CONST blasint M, OPENBLAS_CONST blasint N, OPENBLAS_CONST double alpha, OPENBLAS_CONST double *A, OPENBLAS_CONST blasint lda, double *B, OPENBLAS_CONST blasint ldb){
	cblas_dtrmm(Order, Side, Uplo, TransA, Diag, M, N, alpha, A, lda, B, ldb);
}

template <> 
void trsm(OPENBLAS_CONST enum CBLAS_ORDER Order, OPENBLAS_CONST enum CBLAS_SIDE Side, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,OPENBLAS_CONST enum CBLAS_DIAG Diag, OPENBLAS_CONST blasint M, OPENBLAS_CONST blasint N, OPENBLAS_CONST float alpha, OPENBLAS_CONST float *A, OPENBLAS_CONST blasint lda, float *B, OPENBLAS_CONST blasint ldb){
	cblas_strsm(Order, Side, Uplo, TransA, Diag, M, N, alpha, A, lda, B, ldb);
}

template <> 
void trsm(OPENBLAS_CONST enum CBLAS_ORDER Order, OPENBLAS_CONST enum CBLAS_SIDE Side, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,OPENBLAS_CONST enum CBLAS_DIAG Diag, OPENBLAS_CONST blasint M, OPENBLAS_CONST blasint N, OPENBLAS_CONST double alpha, OPENBLAS_CONST double *A, OPENBLAS_CONST blasint lda, double *B, OPENBLAS_CONST blasint ldb){
	cblas_dtrsm(Order, Side, Uplo, TransA, Diag, M, N, alpha, A, lda, B, ldb);
}

template <> 
void axpby(OPENBLAS_CONST blasint n, OPENBLAS_CONST float alpha, OPENBLAS_CONST float *x, OPENBLAS_CONST blasint incx,OPENBLAS_CONST float beta, float *y, OPENBLAS_CONST blasint incy){
	cblas_saxpby(n, alpha, x, incx, beta, y, incy);
}

template <> 
void axpby(OPENBLAS_CONST blasint n, OPENBLAS_CONST double alpha, OPENBLAS_CONST double *x, OPENBLAS_CONST blasint incx,OPENBLAS_CONST double beta, double *y, OPENBLAS_CONST blasint incy){
	cblas_daxpby(n, alpha, x, incx, beta, y, incy);
}

template <> 
void omatcopy(OPENBLAS_CONST enum CBLAS_ORDER CORDER, OPENBLAS_CONST enum CBLAS_TRANSPOSE CTRANS, OPENBLAS_CONST blasint crows, OPENBLAS_CONST blasint ccols, OPENBLAS_CONST float calpha, OPENBLAS_CONST float *a,OPENBLAS_CONST blasint clda, float *b, OPENBLAS_CONST blasint cldb){
	cblas_somatcopy(CORDER, CTRANS, crows, ccols, calpha, a, clda, b, cldb);
}

template <> 
void omatcopy(OPENBLAS_CONST enum CBLAS_ORDER CORDER, OPENBLAS_CONST enum CBLAS_TRANSPOSE CTRANS, OPENBLAS_CONST blasint crows, OPENBLAS_CONST blasint ccols, OPENBLAS_CONST double calpha, OPENBLAS_CONST double *a,OPENBLAS_CONST blasint clda, double *b, OPENBLAS_CONST blasint cldb){
	cblas_domatcopy(CORDER, CTRANS, crows, ccols, calpha, a, clda, b, cldb);
}

template <> 
void imatcopy(OPENBLAS_CONST enum CBLAS_ORDER CORDER, OPENBLAS_CONST enum CBLAS_TRANSPOSE CTRANS, OPENBLAS_CONST blasint crows, OPENBLAS_CONST blasint ccols, OPENBLAS_CONST float calpha, float *a,OPENBLAS_CONST blasint clda, OPENBLAS_CONST blasint cldb){
	cblas_simatcopy(CORDER, CTRANS, crows, ccols, calpha, a, clda, cldb);
}

template <> 
void imatcopy(OPENBLAS_CONST enum CBLAS_ORDER CORDER, OPENBLAS_CONST enum CBLAS_TRANSPOSE CTRANS, OPENBLAS_CONST blasint crows, OPENBLAS_CONST blasint ccols, OPENBLAS_CONST double calpha, double *a,OPENBLAS_CONST blasint clda, OPENBLAS_CONST blasint cldb){
	cblas_dimatcopy(CORDER, CTRANS, crows, ccols, calpha, a, clda, cldb);
}

}

#endif
