#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "sys/time.h"
#include "assert.h"

// for the OpenMP versions
#include "omp.h"
// for the SIMD versions
#include <xmmintrin.h>


// -------------------------
// double precision versions
// -------------------------

// old baseline implementation using double precision
void layout_fr_double_baseline( double *position, int nodes, 
                       int *es, int edges, int dummy_dim, 
                       int MAX_IT, double temp );

// version using sse (2way double, 128bit registers)
// 2x speedup over double baseline
void layout_fr_double_simd_only( double *position, int nodes, 
                                 int *es, int edges, int dummy_dim,  
                                 int MAX_IT, double temp );


// -------------------------
// single precision versions
// -------------------------

// float version of baseline 
// ~ 1.6x faster than double prec => new baseline
void layout_fr_interleaved_baseline( float *position, int nodes, 
                                     int *es, int edges, int dummy_dim,
                                     int MAX_IT, float temp );

// version using sse (4way float, 128bit registers)
// 1.5x speedup over float baseline
void layout_fr_interleaved_simd_only( float *position, int nodes,
                                      int *es, int edges, int dummy_dim, 
                                      int MAX_IT, float temp );

// OpenMP version
// using private displacement vectors for each thread 
// to avoid race conditions
void layout_fr_interleaved_omp_only( float *position, int nodes,
                                     int *es, int edges, int dummy_dim, 
                                     int MAX_IT, float temp );

// first attempt to go hybrid (simd && OpenMP)
void layout_fr_interleaved( float *position, int nodes,
                            int *es, int edges, int dummy_dim, 
                            int MAX_IT, float temp );


// ---------------------------------
// strided single precision versions
// ---------------------------------

// - idea: use 1D array for easier vectoriztion afterwards.
// is about 5% faster than float baseline
void layout_fr_baseline( float *position, int nodes, 
                         int *es, int edges, int dummy_dim, 
                         int MAX_IT, float temp );

// simd version based on layout_fr_stride
// this version is about 2.9x faster than float baseline and
// thus way better than layout_fr_simd
void layout_fr_simd_only( float *position, int nodes, 
                          int *es, int edges, int dummy_dim, 
                          int MAX_IT, float temp );

// strided ompenmp version
void layout_fr_omp_only( float *position, int nodes, 
                         int *es, int edges, int dummy_dim, 
                         int MAX_IT, float temp );

// second attempt to go hybrid
void layout_fr( float *position, int nodes, 
                int *es, int edges, int dummy_dim, 
                int MAX_IT, float temp );
