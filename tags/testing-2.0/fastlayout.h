#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "sys/time.h"
#include "assert.h"

// for the OpenMP versions
#include "omp.h"
// for the SIMD versions
#include <pmmintrin.h>

// old baseline implementation using double precision
void layout_fr_double(double *position, int nodes, 
                      int *es, int edges, int dummy_dim, 
                      int MAX_IT, double temp );

// float version of baseline 
// ~ 1.6x faster than double prec => new baseline
void layout_fr( float *position, int nodes, 
                int *es, int edges, int dummy_dim,
                int MAX_IT, float temp );

// - version preparing 4 way float vetorization
// - idea: unrolling iloop by 2 and jloop by 2
void layout_fr_2x2( float *position, int nodes, 
                    int *es, int edges, int dummy_dim, 
                    int MAX_IT, float temp );

// TODO: not working
// reinterpret the passed positions for easier vectorization
void layout_fr_stride( float *position, int nodes, 
                       int *es, int edges, int dummy_dim, 
                       int MAX_IT, float temp );

// reorder attractive and repulsive forces loops
// to get scalar replacement
void layout_fr_loop_reorder( float *position, int nodes,
                             int *es, int edges, int dummy_dim,  
                             int MAX_IT, float temp );
// block out i-loop by 2
void layout_fr_block_iloop2( float *position, int nodes, 
                             int *es, int edges, int dummy_dim, 
                             int MAX_IT, float temp );


// OpenMP versions
// using private displacement vectors for each thread 
// to avoid race conditions
void layout_fr_omp( float *position, int nodes,
                    int *es, int edges, int dummy_dim, 
                    int MAX_IT, float temp );

// first attempt to go hybrid
void layout_fr_omp_simd( float *position, int nodes,
                         int *es, int edges, int dummy_dim, 
                         int MAX_IT, float temp );



// SSE version
// version using sse (2way double, 128bit registers)
// 2x speedup over double baseline
void layout_fr_double_simd( double *position, int nodes, 
                            int *es, int edges, int dummy_dim,  
                            int MAX_IT, double temp );

// version using sse (4way float, 128bit registers)
// 1.5x speedup over float baseline
void layout_fr_simd( float *position, int nodes,
                     int *es, int edges, int dummy_dim, 
                     int MAX_IT, float temp );

// alternative approach to get 4way float sse version with
// better speedup
// TODO: not yet working
void layout_fr_simd_jloop( float *position, int nodes,
                           int *es, int edges, int dummy_dim, 
                           int MAX_IT, float temp );

// working 4way float version
// => unroles the inner most j loop by 4
// 1.6x speedup over float baseline
void layout_fr_simd_1x4( float *position, int nodes,
                         int *es, int edges, int dummy_dim, 
                         int MAX_IT, float temp );
