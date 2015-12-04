/*
EAKMeans is a fast Exact K-means library written in C++ with 
command-line interface, shared library + header files and 
Python bindings

Copyright (c) 2015 Idiap Research Institute, http://www.idiap.ch/
Written by James Newling <jnewling@idiap.ch>

This file is part of EAKMeans.

EAKMeans is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3 as
published by the Free Software Foundation.

EAKMeans is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with EAKMeans. If not, see <http://www.gnu.org/licenses/>.



*/

#ifndef BLASTEMPLATES_H 
#define BLASTEMPLATES_H 
#include <cblas.h> 

namespace wblas{
template <typename TFloatType>
TFloatType dot(OPENBLAS_CONST blasint n, OPENBLAS_CONST TFloatType *x, OPENBLAS_CONST blasint incx, OPENBLAS_CONST TFloatType *y, OPENBLAS_CONST blasint incy);

template <typename TFloatType>
TFloatType asum(OPENBLAS_CONST blasint n, OPENBLAS_CONST TFloatType *x, OPENBLAS_CONST blasint incx);

template <typename TFloatType>
TFloatType nrm2(OPENBLAS_CONST blasint N, OPENBLAS_CONST TFloatType *X, OPENBLAS_CONST blasint incX);

template <typename TFloatType>
CBLAS_INDEX i_amax(OPENBLAS_CONST blasint n, OPENBLAS_CONST TFloatType *x, OPENBLAS_CONST blasint incx);

template <typename TFloatType>
void axpy(OPENBLAS_CONST blasint n, OPENBLAS_CONST TFloatType alpha, OPENBLAS_CONST TFloatType *x, OPENBLAS_CONST blasint incx, TFloatType *y, OPENBLAS_CONST blasint incy);

template <typename TFloatType>
void copy(OPENBLAS_CONST blasint n, OPENBLAS_CONST TFloatType *x, OPENBLAS_CONST blasint incx, TFloatType *y, OPENBLAS_CONST blasint incy);

template <typename TFloatType>
void swap(OPENBLAS_CONST blasint n, TFloatType *x, OPENBLAS_CONST blasint incx, TFloatType *y, OPENBLAS_CONST blasint incy);

template <typename TFloatType>
void rot(OPENBLAS_CONST blasint N, TFloatType *X, OPENBLAS_CONST blasint incX, TFloatType *Y, OPENBLAS_CONST blasint incY, OPENBLAS_CONST TFloatType c, OPENBLAS_CONST TFloatType s);

template <typename TFloatType>
void rotg(TFloatType *a, TFloatType *b, TFloatType *c, TFloatType *s);

template <typename TFloatType>
void rotm(OPENBLAS_CONST blasint N, TFloatType *X, OPENBLAS_CONST blasint incX, TFloatType *Y, OPENBLAS_CONST blasint incY, OPENBLAS_CONST TFloatType *P);

template <typename TFloatType>
void rotmg(TFloatType *d1, TFloatType *d2, TFloatType *b1, OPENBLAS_CONST TFloatType b2, TFloatType *P);

template <typename TFloatType>
void scal(OPENBLAS_CONST blasint N, OPENBLAS_CONST TFloatType alpha, TFloatType *X, OPENBLAS_CONST blasint incX);

template <typename TFloatType>
void gemv(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_TRANSPOSE trans, OPENBLAS_CONST blasint m, OPENBLAS_CONST blasint n, OPENBLAS_CONST TFloatType alpha, OPENBLAS_CONST TFloatType *a, OPENBLAS_CONST blasint lda, OPENBLAS_CONST TFloatType *x, OPENBLAS_CONST blasint incx, OPENBLAS_CONST TFloatType beta, TFloatType *y, OPENBLAS_CONST blasint incy);

template <typename TFloatType>
void ger(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST blasint M, OPENBLAS_CONST blasint N, OPENBLAS_CONST TFloatType alpha, OPENBLAS_CONST TFloatType *X, OPENBLAS_CONST blasint incX, OPENBLAS_CONST TFloatType *Y, OPENBLAS_CONST blasint incY, TFloatType *A, OPENBLAS_CONST blasint lda);

template <typename TFloatType>
void trsv(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA, OPENBLAS_CONST enum CBLAS_DIAG Diag, OPENBLAS_CONST blasint N, OPENBLAS_CONST TFloatType *A, OPENBLAS_CONST blasint lda, TFloatType *X, OPENBLAS_CONST blasint incX);

template <typename TFloatType>
void trmv(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA, OPENBLAS_CONST enum CBLAS_DIAG Diag, OPENBLAS_CONST blasint N, OPENBLAS_CONST TFloatType *A, OPENBLAS_CONST blasint lda, TFloatType *X, OPENBLAS_CONST blasint incX);

template <typename TFloatType>
void syr(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST blasint N, OPENBLAS_CONST TFloatType alpha, OPENBLAS_CONST TFloatType *X, OPENBLAS_CONST blasint incX, TFloatType *A, OPENBLAS_CONST blasint lda);

template <typename TFloatType>
void syr2(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST blasint N, OPENBLAS_CONST TFloatType alpha, OPENBLAS_CONST TFloatType *X, OPENBLAS_CONST blasint incX, OPENBLAS_CONST TFloatType *Y, OPENBLAS_CONST blasint incY, TFloatType *A, OPENBLAS_CONST blasint lda);

template <typename TFloatType>
void gbmv(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA, OPENBLAS_CONST blasint M, OPENBLAS_CONST blasint N, OPENBLAS_CONST blasint KL, OPENBLAS_CONST blasint KU, OPENBLAS_CONST TFloatType alpha, OPENBLAS_CONST TFloatType *A, OPENBLAS_CONST blasint lda, OPENBLAS_CONST TFloatType *X, OPENBLAS_CONST blasint incX, OPENBLAS_CONST TFloatType beta, TFloatType *Y, OPENBLAS_CONST blasint incY);

template <typename TFloatType>
void sbmv(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST blasint N, OPENBLAS_CONST blasint K, OPENBLAS_CONST TFloatType alpha, OPENBLAS_CONST TFloatType *A, OPENBLAS_CONST blasint lda, OPENBLAS_CONST TFloatType *X, OPENBLAS_CONST blasint incX, OPENBLAS_CONST TFloatType beta, TFloatType *Y, OPENBLAS_CONST blasint incY);

template <typename TFloatType>
void tbmv(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA, OPENBLAS_CONST enum CBLAS_DIAG Diag, OPENBLAS_CONST blasint N, OPENBLAS_CONST blasint K, OPENBLAS_CONST TFloatType *A, OPENBLAS_CONST blasint lda, TFloatType *X, OPENBLAS_CONST blasint incX);

template <typename TFloatType>
void tbsv(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA, OPENBLAS_CONST enum CBLAS_DIAG Diag, OPENBLAS_CONST blasint N, OPENBLAS_CONST blasint K, OPENBLAS_CONST TFloatType *A, OPENBLAS_CONST blasint lda, TFloatType *X, OPENBLAS_CONST blasint incX);

template <typename TFloatType>
void tpmv(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA, OPENBLAS_CONST enum CBLAS_DIAG Diag, OPENBLAS_CONST blasint N, OPENBLAS_CONST TFloatType *Ap, TFloatType *X, OPENBLAS_CONST blasint incX);

template <typename TFloatType>
void tpsv(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA, OPENBLAS_CONST enum CBLAS_DIAG Diag, OPENBLAS_CONST blasint N, OPENBLAS_CONST TFloatType *Ap, TFloatType *X, OPENBLAS_CONST blasint incX);

template <typename TFloatType>
void symv(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST blasint N, OPENBLAS_CONST TFloatType alpha, OPENBLAS_CONST TFloatType *A, OPENBLAS_CONST blasint lda, OPENBLAS_CONST TFloatType *X, OPENBLAS_CONST blasint incX, OPENBLAS_CONST TFloatType beta, TFloatType *Y, OPENBLAS_CONST blasint incY);

template <typename TFloatType>
void spmv(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST blasint N, OPENBLAS_CONST TFloatType alpha, OPENBLAS_CONST TFloatType *Ap, OPENBLAS_CONST TFloatType *X, OPENBLAS_CONST blasint incX, OPENBLAS_CONST TFloatType beta, TFloatType *Y, OPENBLAS_CONST blasint incY);

template <typename TFloatType>
void spr(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST blasint N, OPENBLAS_CONST TFloatType alpha, OPENBLAS_CONST TFloatType *X, OPENBLAS_CONST blasint incX, TFloatType *Ap);

template <typename TFloatType>
void spr2(OPENBLAS_CONST enum CBLAS_ORDER order, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST blasint N, OPENBLAS_CONST TFloatType alpha, OPENBLAS_CONST TFloatType *X, OPENBLAS_CONST blasint incX, OPENBLAS_CONST TFloatType *Y, OPENBLAS_CONST blasint incY, TFloatType *A);

template <typename TFloatType>
void gemm(OPENBLAS_CONST enum CBLAS_ORDER Order, OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA, OPENBLAS_CONST enum CBLAS_TRANSPOSE TransB, OPENBLAS_CONST blasint M, OPENBLAS_CONST blasint N, OPENBLAS_CONST blasint K, OPENBLAS_CONST TFloatType alpha, OPENBLAS_CONST TFloatType *A, OPENBLAS_CONST blasint lda, OPENBLAS_CONST TFloatType *B, OPENBLAS_CONST blasint ldb, OPENBLAS_CONST TFloatType beta, TFloatType *C, OPENBLAS_CONST blasint ldc);

template <typename TFloatType>
void symm(OPENBLAS_CONST enum CBLAS_ORDER Order, OPENBLAS_CONST enum CBLAS_SIDE Side, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST blasint M, OPENBLAS_CONST blasint N, OPENBLAS_CONST TFloatType alpha, OPENBLAS_CONST TFloatType *A, OPENBLAS_CONST blasint lda, OPENBLAS_CONST TFloatType *B, OPENBLAS_CONST blasint ldb, OPENBLAS_CONST TFloatType beta, TFloatType *C, OPENBLAS_CONST blasint ldc);

template <typename TFloatType>
void syrk(OPENBLAS_CONST enum CBLAS_ORDER Order, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST enum CBLAS_TRANSPOSE Trans, OPENBLAS_CONST blasint N, OPENBLAS_CONST blasint K, OPENBLAS_CONST TFloatType alpha, OPENBLAS_CONST TFloatType *A, OPENBLAS_CONST blasint lda, OPENBLAS_CONST TFloatType beta, TFloatType *C, OPENBLAS_CONST blasint ldc);

template <typename TFloatType>
void syr2k(OPENBLAS_CONST enum CBLAS_ORDER Order, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST enum CBLAS_TRANSPOSE Trans, OPENBLAS_CONST blasint N, OPENBLAS_CONST blasint K, OPENBLAS_CONST TFloatType alpha, OPENBLAS_CONST TFloatType *A, OPENBLAS_CONST blasint lda, OPENBLAS_CONST TFloatType *B, OPENBLAS_CONST blasint ldb, OPENBLAS_CONST TFloatType beta, TFloatType *C, OPENBLAS_CONST blasint ldc);

template <typename TFloatType>
void trmm(OPENBLAS_CONST enum CBLAS_ORDER Order, OPENBLAS_CONST enum CBLAS_SIDE Side, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA, OPENBLAS_CONST enum CBLAS_DIAG Diag, OPENBLAS_CONST blasint M, OPENBLAS_CONST blasint N, OPENBLAS_CONST TFloatType alpha, OPENBLAS_CONST TFloatType *A, OPENBLAS_CONST blasint lda, TFloatType *B, OPENBLAS_CONST blasint ldb);

template <typename TFloatType>
void trsm(OPENBLAS_CONST enum CBLAS_ORDER Order, OPENBLAS_CONST enum CBLAS_SIDE Side, OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA, OPENBLAS_CONST enum CBLAS_DIAG Diag, OPENBLAS_CONST blasint M, OPENBLAS_CONST blasint N, OPENBLAS_CONST TFloatType alpha, OPENBLAS_CONST TFloatType *A, OPENBLAS_CONST blasint lda, TFloatType *B, OPENBLAS_CONST blasint ldb);

template <typename TFloatType>
void axpby(OPENBLAS_CONST blasint n, OPENBLAS_CONST TFloatType alpha, OPENBLAS_CONST TFloatType *x, OPENBLAS_CONST blasint incx, OPENBLAS_CONST TFloatType beta, TFloatType *y, OPENBLAS_CONST blasint incy);

template <typename TFloatType>
void omatcopy(OPENBLAS_CONST enum CBLAS_ORDER CORDER, OPENBLAS_CONST enum CBLAS_TRANSPOSE CTRANS, OPENBLAS_CONST blasint crows, OPENBLAS_CONST blasint ccols, OPENBLAS_CONST TFloatType calpha, OPENBLAS_CONST TFloatType *a, OPENBLAS_CONST blasint clda, TFloatType *b, OPENBLAS_CONST blasint cldb);

template <typename TFloatType>
void imatcopy(OPENBLAS_CONST enum CBLAS_ORDER CORDER, OPENBLAS_CONST enum CBLAS_TRANSPOSE CTRANS, OPENBLAS_CONST blasint crows, OPENBLAS_CONST blasint ccols, OPENBLAS_CONST TFloatType calpha, TFloatType *a, OPENBLAS_CONST blasint clda, OPENBLAS_CONST blasint cldb);

} 

#endif


