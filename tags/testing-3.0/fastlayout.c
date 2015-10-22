#include "fastlayout.h"

// aligned memory allocation
#include "memory.c"

/* helper function to get time duration from start and end values
 */
inline double gettime( struct timeval start, struct timeval end ) {
  return (double)((end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)/1e6);
}

/* return a random number in [0, max)
 */
inline float get_random_positions( float max ) {
  return max * (rand() / (float) ((unsigned)RAND_MAX+1)) - max/2;
}

/* helper function to randomly initialize positions
 */
void initialize( float* positions, int num_nodes, int dimension ) {
  float sqrt_nodes = sqrt(num_nodes);
  unsigned i;
  for( i = 0; i < num_nodes*dimension; ++i ) {
    positions[i] = get_random_positions(sqrt_nodes);
  }
}

//////////////////////////////////////////
// Layouts in double precision
// ---------------------------


/* layout fruchterman reingold
* input: 
* - problem size: #nodes * dimension
* - edge list
*   2d array, first dimension source, second dim dest
* - MAX_IT:  number of iterations
* - temp:    some sort of temperature 
*/
// double precision baseline
void layout_fr_double( double *position, int nodes, 
                       int *es, int edges, int dummy_dim, 
                       int MAX_IT, double temp ) {
  int n, i, j;
  int source, target;
  
  double dx, dy, len;
  double *dplx = (double*) aligned_malloc(nodes*sizeof(double));
  double *dply = (double*) aligned_malloc(nodes*sizeof(double));
  
  double real_dx, real_dy;
  double difftemp = temp / MAX_IT;
  
  // first, set random positions
  const unsigned int dim = 2;
  nodes /= dim;
  double sqrt_nodes = sqrt(nodes);
  for( i = 0; i < nodes*dim; ++i ) {
    position[i] = sqrt_nodes * (rand() / (double) ((unsigned)RAND_MAX+1)) - sqrt_nodes/2;
  }
  
  // do this for some iterations
  for( n=0; n < MAX_IT; ++n ) {
    
    // clear displacement vectors
    for( i=0; i < nodes; ++i ) {
      dplx[i] = 0.;
      dply[i] = 0.;
    }
    
    // repulsive forces
    for( i = 0; i < nodes; ++i ) {
      for( j = i+1; j < nodes; ++j ) {
        dx = position[i*dim] - position[j*dim];
        dy = position[i*dim+1] - position[j*dim+1];
        len = dx*dx + dy*dy;
        
        // avoid division by zero
        if( len == 0. ) {
          dx = rand() * 1e-9;
          dy = rand() * 1e-9;
          len = dx*dx + dy*dy;
        }
        
        dplx[i] += dx/len;
        dply[i] += dy/len;
        dplx[j] -= dx/len;
        dply[j] -= dy/len;
      }
    }
    
    // attractive forces
    for( i = 0; i < edges; ++i ) {
      source = es[i*dim];
      target = es[i*dim+1];
      
      dx = position[source*dim] - position[target*dim];
      dy = position[source*dim+1] - position[target*dim+1];
      len = sqrt(dx*dx + dy*dy);
      
      dplx[source] -= dx * len;
      dply[source] -= dy * len;
      dplx[target] += dx * len;
      dply[target] += dy * len;
    }
    
    // finally update the positions
    for( i=0; i < nodes; ++i ) {
      dx = dplx[i] + rand() * 1e-9;
      dy = dply[i] + rand() * 1e-9;
      len = sqrt( dx*dx + dy*dy );
      
      real_dx = ( fabs(dx) < temp ) ? dx : temp;
      real_dy = ( fabs(dy) < temp ) ? dy : temp;
      
      // avoid division by zero
      if( len > 0. ) {
        position[i*dim] += (dx / len) * real_dx;
        position[i*dim+1] += (dy / len) * real_dy;
      }
    }
    
    temp -= difftemp;
  } // end of iteration loop
  
  // free allocated memory
  aligned_free( dplx );
  aligned_free( dply );
}

/* layout fruchterman reingold
* input: 
* - problem size: #nodes * dimension
* - edge list
*   2d array, first dimension source, second dim dest
* - MAX_IT:  number of iterations
* - temp:    some sort of temperature 
*/
// working 2-way double precision SSE version
// 2x speedup over double precision baseline
void layout_fr_double_simd( double *position, int nodes, 
                            int *es, int edges, int dummy_dim,  
                            int MAX_IT, double temp ) {
  int n, i, j;
  int source, target;
  int STRIDE = 2;
  
  double dx, dy, dx2, dy2, len, len2;
  double dxi, dyi, pxi, pyi;
  double *dplx = (double*) aligned_malloc(nodes*sizeof(double));
  double *dply = (double*) aligned_malloc(nodes*sizeof(double));
  
  double real_dx, real_dy;
  double difftemp = temp / MAX_IT;
  
  __m128d di, dii, pi, pii, pj;
  __m128d d, d2, l, l2;
  
  double *tmp = (double*) aligned_malloc(4*sizeof(double));
  
  // first, set random positions
  const unsigned int dim = 2;
  nodes /= dim;
  double sqrt_nodes = sqrt(nodes);
  for( i = 0; i < nodes*dim; ++i ) {
    position[i] = sqrt_nodes * (rand() / (double) ((unsigned)RAND_MAX+1)) - sqrt_nodes/2;
  }
  
  // do this for some iterations
  for( n=0; n < MAX_IT; ++n ) {
    // repulsive forces
    for( i = 0; (i+STRIDE-1) < nodes; i+=STRIDE ) {
      // clear displacement vectors
      di = _mm_setzero_pd();
      dii = _mm_setzero_pd();
      
      // load positions of nodes i and i+1 
      pi = _mm_loadu_pd(&position[i*dim]);
      pii = _mm_loadu_pd(&position[(i+1)*dim]);
      
      // special treatment for the fist j
      // this needs only be processed with pi and not pii
      j = i+1;
      if( j < nodes ) {
        pj = _mm_loadu_pd( &position[j*dim] );
        
        d = _mm_sub_pd( pi, pj );
        _mm_store_pd( tmp, _mm_mul_pd(d, d) );
        len = tmp[0] + tmp[1];
        
        // avoid division by zero
        if( len == 0. ) {
          dx = rand() * 1e-9;
          dy = rand() * 1e-9;
          len = dx*dx + dy*dy;
        }
        
        l = _mm_load_pd1( &len );
        d = _mm_div_pd( d, l );
        di = _mm_add_pd( di, d );
        
        _mm_store_pd(tmp, d);
        dplx[j] -= tmp[0];
        dply[j] -= tmp[1];
      }
      // do the other j's
      for(j = j+2 ; j < nodes; ++j ) {
        // load position of node j
        pj = _mm_loadu_pd(&position[j*dim]);
        
        // distance to nodes i and i+1
        d = _mm_sub_pd( pi, pj );
        d2 = _mm_sub_pd(pii, pj);
        // lenght of distance
        _mm_store_pd(tmp, _mm_mul_pd(d, d));
        _mm_store_pd(tmp + 2, _mm_mul_pd(d2, d2));
        len  = tmp[0] + tmp[1];
        len2 = tmp[2] + tmp[3];
        
        // avoid division by zero
        if( len == 0. ) {
          dx = rand() * 1e-9;
          dy = rand() * 1e-9;
          len = dx*dx + dy*dy;
        }
        if( len2 == 0. ) {
          dx2 = rand() * 1e-9;
          dy2 = rand() * 1e-9;
          len2 = dx2*dx2 + dy2*dy2;
        }
        
        l = _mm_load_pd1( &len );
        l2 = _mm_load_pd1( &len2 );
        
        d = _mm_div_pd( d, l);
        d2 = _mm_div_pd( d2, l2 );
        // update displacement vector for i and i+1
        di = _mm_add_pd( di, d );
        dii = _mm_add_pd( dii, d2 );
        
        d = _mm_add_pd(d, d2);
        _mm_store_pd(tmp, d);
        // update displacement vector for j
        dplx[j] -= tmp[0];
        dply[j] -= tmp[1];
      }
      _mm_store_pd( tmp, di );
      _mm_store_pd( tmp + 2, dii );
      dplx[i] = tmp[0];
      dplx[i+1] = tmp[2];
      dply[i] = tmp[1];
      dply[i+1] = tmp[3];
    }
    for( ; i < nodes; ++i ) {
      // clear displacement vectors
      dxi = 0.;
      dyi = 0.;
      
      pxi = position[i*dim];
      pyi = position[i*dim+1];
      
      for( j = i+1; j < nodes; ++j ) {
        dx = pxi - position[j*dim];
        dy = pyi - position[j*dim+1];
        len = dx*dx + dy*dy;
        
        // avoid division by zero
        if( len == 0. ) {
          dx = rand() * 1e-9;
          dy = rand() * 1e-9;
          len = dx*dx + dy*dy;
        }
        
        dxi += dx/len;
        dyi += dy/len;
        dplx[j] -= dx/len;
        dply[j] -= dy/len;   
      }
      dplx[i] = dxi;
      dply[i] = dyi;
    }
    
    // attractive forces
    for( i = 0; i < edges; ++i ) {
      source = es[i*dim];
      target = es[i*dim+1];
      
      dx = position[source*dim] - position[target*dim];
      dy = position[source*dim+1] - position[target*dim+1];
      len = sqrt(dx*dx + dy*dy);
      
      dplx[source] -= dx * len;
      dply[source] -= dy * len;
      dplx[target] += dx * len;
      dply[target] += dy * len;
    }
    
    // finally update the positions
    for( i=0; i < nodes; ++i ) {
      dx = dplx[i] + rand() * 1e-9;
      dy = dply[i] + rand() * 1e-9;
      len = sqrt( dx*dx + dy*dy );

      real_dx = ( fabs(dx) < temp ) ? dx : temp;
      real_dy = ( fabs(dy) < temp ) ? dy : temp;
      
      // avoid division by zero
      if( len > 0. ) {
        position[i*dim] += (dx / len) * real_dx;
        position[i*dim+1] += (dy / len) * real_dy;
      }
    }
    temp -= difftemp;
  } // end of iteration loop

  
  // free allocated memory
  aligned_free( dplx );
  aligned_free( dply );
  aligned_free( tmp );
}

////////////////////////////////////////////////////
// Layouts in single precision
// ---------------------------

/* layout fruchterman reingold
* input: 
* - problem size: #nodes * dimension 
* - edge list
*   2d array, first dimension source, second dim dest
* - MAX_IT:  number of iterations
* - temp:    some sort of temperature 
*/
// float version of 'layout_fr'
void layout_fr( float *position, int nodes, 
                int *es, int edges, int dummy_dim,
                int MAX_IT, float temp ) {
  int n, i, j;
  int source, target;
  
  float dx, dy, len;
  float *dplx = (float*) aligned_malloc(nodes*sizeof(float));
  float *dply = (float*) aligned_malloc(nodes*sizeof(float));
  
  float real_dx, real_dy;
  float difftemp = temp / MAX_IT;
  
  // first, set random positions
  const unsigned int dim = 2;
  nodes /= dim;
  initialize( position, nodes, dim );
  
  // do this for some iterations
  for( n=0; n < MAX_IT; ++n ) {
    
    // clear displacement vectors
    for( i=0; i < nodes; ++i ) {
      dplx[i] = 0.;
      dply[i] = 0.;
    }
    
    // repulsive forces
    for( i = 0; i < nodes; ++i ) {
      for( j = i+1; j < nodes; ++j ) {
        dx = position[i*dim] - position[j*dim];
        dy = position[i*dim+1] - position[j*dim+1];
        len = dx*dx + dy*dy;
        
        // avoid division by zero
        if( len == 0. ) {
          dx = rand() * 1e-9;
          dy = rand() * 1e-9;
          len = dx*dx + dy*dy;
        }
        
        dplx[i] += dx/len;
        dply[i] += dy/len;
        dplx[j] -= dx/len;
        dply[j] -= dy/len;
      }
    }
    
    // attractive forces
    for( i = 0; i < edges; ++i ) {
      source = es[i*dim];
      target = es[i*dim+1];
      
      dx = position[source*dim] - position[target*dim];
      dy = position[source*dim+1] - position[target*dim+1];
      len = sqrt(dx*dx + dy*dy);
      
      dplx[source] -= dx * len;
      dply[source] -= dy * len;
      dplx[target] += dx * len;
      dply[target] += dy * len;
    }
    
    // finally update the positions
    for( i=0; i < nodes; ++i ) {
      dx = dplx[i] + rand() * 1e-9;
      dy = dply[i] + rand() * 1e-9;
      len = sqrt( dx*dx + dy*dy );
      
      real_dx = ( fabs(dx) < temp ) ? dx : temp;
      real_dy = ( fabs(dy) < temp ) ? dy : temp;
      
      // avoid division by zero
      if( len > 0. ) {
        position[i*dim] += (dx / len) * real_dx;
        position[i*dim+1] += (dy / len) * real_dy;
      }
    }
    
    temp -= difftemp;
  } // end of iteration loop
  
  // free allocated memory
  aligned_free( dplx );
  aligned_free( dply );
}

/* layout fruchterman reingold
* input: 
* - problem size: #nodes * dimension 
* - edge list
*   2d array, first dimension source, second dim dest
* - MAX_IT:  number of iterations
* - temp:    some sort of temperature 
*/
// working single precision 4way float version
// 1.5x speedup compared to float baseversion
void layout_fr_simd( float *position, int nodes,
                     int *es, int edges, int dummy_dim, 
                     int MAX_IT, float temp ) {
  // NOTE
  // this version seems broken 
  // it gives very different pictures for the test graph
  int n, i, j;
  int source, target;
  int STRIDE = 2;
  
  float dx, dy, dx2, dy2, len, len2;
  float dxi, dyi, pxi, pyi;
  float *dplx = (float*) aligned_malloc(nodes*sizeof(float));
  float *dply = (float*) aligned_malloc(nodes*sizeof(float));
  
  float real_dx, real_dy;
  float difftemp = temp / MAX_IT;
  
  __m128 di, pi, pj;
  __m128 d, l;
  float *tmp = (float*) aligned_malloc(4*sizeof(float));
  
  // first, set random positions
  const unsigned int dim = 2;
  nodes /= dim;
  initialize( position, nodes, dim );
  
  // do this for some iterations
  for( n=0; n < MAX_IT; ++n ) {
    
    // clear displacement vectors
    for( i = 0; i < nodes; ++i ) {
      dplx[i] = 0.;
      dply[i] = 0.;
    }

    // repulsive forces
    for( i = 0; (i+STRIDE-1) < nodes; i+=STRIDE ) {
      // clear displacement vectors
      di = _mm_setzero_ps();
      
      // load positions of nodes i and i+1 
      pi = _mm_loadu_ps(&position[i*dim]);
      pj = pi;
      
      // special treatment for the fist j
      // this needs only be processed with pi and not pii
      j = i+1;
      if( j < nodes ) {
        // load p^j in the lower two slots
        pj = _mm_loadl_pi(pj, (__m64*)&position[j*dim] );
        
        // d = [p^i_x - p^j_x, p^i_y p^j_y, 0, 0]
        d = _mm_sub_ps( pi, pj );
        _mm_store_ps( tmp, _mm_mul_ps(d, d) );
        assert( tmp[2] == 0 );
        assert( tmp[3] == 0 );
        len = tmp[0] + tmp[1];
        
        // avoid division by zero
        if( len == 0. ) {
          dx = rand() * 1e-9;
          dy = rand() * 1e-9;
          len = dx*dx + dy*dy;
        }
        
        l = _mm_set1_ps( len );
        d = _mm_div_ps( d, l );
        di = _mm_add_ps( di, d );
        
        _mm_store_ps(tmp, d);
        assert( tmp[2] == 0 );
        assert( tmp[3] == 0 );
        
        // NOTE
        // as j is always bigger than i, dplx/y[j] has always
        // been cleared by _mm_setzero_ps 
        dplx[j] -= tmp[0];
        dply[j] -= tmp[1];
      }
      // do the other j's
      for( j = i+2; j < nodes; ++j ) {
        // load position of node j
        pj = _mm_loadl_pi(pj, (__m64*)&position[j*dim]);
        pj = _mm_loadh_pi(pj, (__m64*)&position[j*dim]);
        
        // distance to nodes i and i+1
        d = _mm_sub_ps( pi, pj );
        
        // lenght of distance
        _mm_store_ps(tmp, _mm_mul_ps(d, d));
        len  = tmp[0] + tmp[1];
        len2 = tmp[2] + tmp[3];
        
        // avoid division by zero
        if( len == 0. ) {
          dx = rand() * 1e-9;
          dy = rand() * 1e-9;
          len = dx*dx + dy*dy;
        }
        if( len2 == 0. ) {
          dx2 = rand() * 1e-9;
          dy2 = rand() * 1e-9;
          len2 = dx2*dx2 + dy2*dy2;
        }
        
        //     LSB
        //     |
        // l = [len, len, len2, len2]
        l = _mm_set_ps(len2, len2, len, len);
        d = _mm_div_ps( d, l);
        
        // update displacement vector for i and i+1
        di = _mm_add_ps( di, d );
        
        _mm_store_ps(tmp, d);
        // update displacement vector for j
        dplx[j] -= (tmp[0] + tmp[2]);
        dply[j] -= (tmp[1] + tmp[3]);
      }
      _mm_store_ps( tmp, di );
      dplx[i]   += tmp[0];
      dplx[i+1] += tmp[2];
      dply[i]   += tmp[1];
      dply[i+1] += tmp[3];
    }
    for( ; i < nodes; ++i ) {
      // clear displacement vectors
      dxi = 0.;
      dyi = 0.;
      
      pxi = position[i*dim];
      pyi = position[i*dim+1];
      
      for( j = i+1; j < nodes; ++j ) {
        dx = pxi - position[j*dim];
        dy = pyi - position[j*dim+1];
        len = dx*dx + dy*dy;
        
        // avoid division by zero
        if( len == 0. ) {
          dx = rand() * 1e-9;
          dy = rand() * 1e-9;
          len = dx*dx + dy*dy;
        }
        
        dxi += dx/len;
        dyi += dy/len;
        dplx[j] -= dx/len;
        dply[j] -= dy/len;   
      }
      dplx[i] += dxi;
      dply[i] += dyi;
    }
    
    // attractive forces
    for( i = 0; i < edges; ++i ) {
      source = es[i*dim];
      target = es[i*dim+1];
      
      dx = position[source*dim] - position[target*dim];
      dy = position[source*dim+1] - position[target*dim+1];
      len = sqrt(dx*dx + dy*dy);
      
      dplx[source] -= dx * len;
      dply[source] -= dy * len;
      dplx[target] += dx * len;
      dply[target] += dy * len;
    }
    
    // finally update the positions
    for( i=0; i < nodes; ++i ) {
      dx = dplx[i] + rand() * 1e-9;
      dy = dply[i] + rand() * 1e-9;
      len = sqrt( dx*dx + dy*dy );

      real_dx = ( fabs(dx) < temp ) ? dx : temp;
      real_dy = ( fabs(dy) < temp ) ? dy : temp;
      
      // avoid division by zero
      if( len > 0. ) {
        position[i*dim] += (dx / len) * real_dx;
        position[i*dim+1] += (dy / len) * real_dy;
      }
    }
    temp -= difftemp;
  } // end of iteration loop

  
  // free allocated memory
  aligned_free( dplx );
  aligned_free( dply );
  aligned_free( tmp );
}

/* layout fruchterman reingold
* input: 
* - problem size: #nodes * dimension
* - edge list
*   2d array, first dimension source, second dim dest
* - MAX_IT:  number of iterations
* - temp:    some sort of temperature 
*/
void layout_fr_omp( float *position, int nodes,
                    int *es, int edges, int dummy_dim, 
                    int MAX_IT, float temp ) {
  struct timeval start, end;
  gettimeofday(&start, NULL);
  int n, i, j;
  int source, target;
  
  float dx, dy, len;
  float dxi, dyi, real_dx, real_dy, px, py;
  float *dplx = (float*) malloc(nodes*sizeof(float));
  float *dply = (float*) malloc(nodes*sizeof(float));
  
  float *space = (float*) malloc(2*nodes*sizeof(float));
  
  float difftemp = temp / MAX_IT;
  
  // first, set random positions
  const unsigned int dim = 2;
  nodes /= dim;
  initialize( position, nodes, dim );
  
  // do this for some iterations
  for( n=0; n < MAX_IT; ++n ) {
    
    for( i = 0; i < nodes; ++i ) {
      dplx[i] = 0.;
      dply[i] = 0.;
    }
    
    // attractive forces
    for( i = 0; i < edges; ++i ) {
      source = es[i*dim];
      target = es[i*dim+1];
      
      dx = position[source*dim] - position[target*dim];
      dy = position[source*dim+1] - position[target*dim+1];
      len = sqrt(dx*dx + dy*dy);
      
      dplx[source] -= dx * len;
      dply[source] -= dy * len;
      dplx[target] += dx * len;
      dply[target] += dy * len;
    }
    
    // repulsive forces
    #pragma omp parallel shared(position, dplx, dply, temp, nodes) private(i,j, dx, dy, len, dxi, dyi, px, py)
    { 
      float *dx_local = space;
      float *dy_local = space + nodes;
      for(i = 0; i < nodes; ++i ) {
        dx_local[i] = 0.;
        dy_local[i] = 0.;
      }
      
      // NOTE:
      // dynamic scheduling is important!!
      // static scheduling fails because the work in jloop
      // is not equal for all i
      #pragma omp for schedule(dynamic)
      for( i = 0; i < nodes; ++i ) {
        dxi = 0.;
        dyi = 0.;
        
        px = position[i*dim];
        py = position[i*dim+1];
        
        for( j = i+1; j < nodes; ++j ) {
          dx = px - position[j*dim];
          dy = py - position[j*dim+1];
          len = dx*dx + dy*dy;
          
          // avoid division by zero
          if( len == 0. ) {
            dx = rand() * 1e-9;
            dy = rand() * 1e-9;
            len = dx*dx + dy*dy;
          }
          
          dxi += dx/len;
          dyi += dy/len;
          dx_local[j] -= dx/len;
          dy_local[j] -= dy/len;
        }
        dx_local[i] += dxi;
        dy_local[i] += dyi;
      } // end of double for loop
      
      #pragma omp for
      for(i = 0; i < nodes; ++i ) {

        dxi = dplx[i];
        dxi += dx_local[i];
        dyi = dply[i];
        dyi += dy_local[i];

        dxi += rand() * 1e-9;
        dyi += rand() * 1e-9;
        len = sqrt( dxi*dxi + dyi*dyi );
        
        real_dx = ( fabs(dxi) < temp ) ? dxi : temp;
        real_dy = ( fabs(dyi) < temp ) ? dyi : temp;
        
        // avoid division by zero
        if( len > 0. ) {
          position[i*dim] += (dxi / len) * real_dx;
          position[i*dim+1] += (dyi / len) * real_dy;
        }
      }
    }
    temp -= difftemp;
  } // end of iteration loop
  
  // free allocated memory
  free( dplx );
  free( dply );
  free( space );
  
  gettimeofday(&end, NULL);
  double time_elapsed = gettime(start,end);
  printf("Layouting this graph took %1.4f sec\n", time_elapsed);
}

/* layout fruchterman reingold
* input: 
* - problem size: #nodes * dimension 
* - edge list
*   2d array, first dimension source, second dim dest
* - MAX_IT:  number of iterations
* - temp:    some sort of temperature 
*/
// first attempt to go hybrid
void layout_fr_omp_simd( float *position, int nodes,
                         int *es, int edges, int dummy_dim, 
                         int MAX_IT, float temp ) {
  struct timeval start, end;
  gettimeofday(&start, NULL);
  int n, i, j;
  int source, target;
  
  float dx, dy, len;
  float dxi, dyi, real_dx, real_dy, px, py;
  float *dplx = (float*) malloc(nodes*sizeof(float));
  float *dply = (float*) malloc(nodes*sizeof(float));
  float *space = (float*) malloc(2*nodes*sizeof(float));
  
  float difftemp = temp / MAX_IT;
  
  const unsigned int STRIDE = 2;
  float dx2, dy2,len2, pxi, pyi;
  __m128 di, pi, pj;
  __m128 d, l;
  
  // first, set random positions
  const unsigned int dim = 2;
  nodes /= dim;
  initialize( position, nodes, dim );
  
  // do this for some iterations
  for( n=0; n < MAX_IT; ++n ) {
    
    for( i = 0; i < nodes; ++i ) {
      dplx[i] = 0.;
      dply[i] = 0.;
    }
    
    // attractive forces
    for( i = 0; i < edges; ++i ) {
      source = es[i*dim];
      target = es[i*dim+1];
      
      dx = position[source*dim] - position[target*dim];
      dy = position[source*dim+1] - position[target*dim+1];
      len = sqrt(dx*dx + dy*dy);
      
      dplx[source] -= dx * len;
      dply[source] -= dy * len;
      dplx[target] += dx * len;
      dply[target] += dy * len;
    }
    
    // repulsive forces
    #pragma omp parallel shared(position, dplx, dply, temp, nodes) private(i,j, dx, dy, len, len2, dx2, dy2, pxi, pyi, dxi, dyi, px, py, di, pi, pj, d, l)
    { 
      float *dx_local = space;
      float *dy_local = space + nodes;
      for(i = 0; i < nodes; ++i ) {
        dx_local[i] = 0.;
        dy_local[i] = 0.;
      }
      float *tmp = (float*) aligned_malloc(4*sizeof(float));
      
      // NOTE:
      // dynamic scheduling is important!!
      // static scheduling fails because the work in jloop
      // is not equal for all i
      #pragma omp for schedule(dynamic)
      for( i = 0; i < (nodes-STRIDE+1); i+=STRIDE ) {
        // clear displacement vectors
        di = _mm_setzero_ps();
        
        // load positions of nodes i and i+1 
        pi = _mm_loadu_ps(&position[i*dim]);
        pj = pi;
        
        // special treatment for the fist j
        // this needs only be processed with pi and not pii
        j = i+1;
        if( j < nodes ) {
          // load p^j in the lower two slots
          pj = _mm_loadl_pi(pj, (__m64*)&position[j*dim] );
          
          // d = [p^i_x - p^j_x, p^i_y p^j_y, 0, 0]
          d = _mm_sub_ps( pi, pj );
          _mm_store_ps( tmp, _mm_mul_ps(d, d) );
          assert( tmp[2] == 0 );
          assert( tmp[3] == 0 );
          len = tmp[0] + tmp[1];
          
          // avoid division by zero
          if( len == 0. ) {
            dx = rand() * 1e-9;
            dy = rand() * 1e-9;
            len = dx*dx + dy*dy;
          }
          
          l = _mm_set1_ps( len );
          d = _mm_div_ps( d, l );
          di = _mm_add_ps( di, d );
          
          _mm_store_ps(tmp, d);
          assert( tmp[2] == 0 );
          assert( tmp[3] == 0 );
          
          // NOTE
          // as j is always bigger than i, dplx/y[j] has always
          // been cleared by _mm_setzero_ps 
          dx_local[j] -= tmp[0];
          dy_local[j] -= tmp[1];
        }
        // do the other j's
        for( j = i+2; j < nodes; ++j ) {
          // load position of node j
          pj = _mm_loadl_pi(pj, (__m64*)&position[j*dim]);
          pj = _mm_loadh_pi(pj, (__m64*)&position[j*dim]);
          
          // distance to nodes i and i+1
          d = _mm_sub_ps( pi, pj );
          
          // lenght of distance
          _mm_store_ps(tmp, _mm_mul_ps(d, d));
          len  = tmp[0] + tmp[1];
          len2 = tmp[2] + tmp[3];
          
          // avoid division by zero
          if( len == 0. ) {
            dx = rand() * 1e-9;
            dy = rand() * 1e-9;
            len = dx*dx + dy*dy;
          }
          if( len2 == 0. ) {
            dx2 = rand() * 1e-9;
            dy2 = rand() * 1e-9;
            len2 = dx2*dx2 + dy2*dy2;
          }
          
          //     LSB
          //     |
          // l = [len, len, len2, len2]
          l = _mm_set_ps(len2, len2, len, len);
          d = _mm_div_ps( d, l);
          
          // update displacement vector for i and i+1
          di = _mm_add_ps( di, d );
          
          _mm_store_ps(tmp, d);
          // update displacement vector for j
          dx_local[j] -= (tmp[0] + tmp[2]);
          dy_local[j] -= (tmp[1] + tmp[3]);
        }
        _mm_store_ps( tmp, di );
        dx_local[i]   += tmp[0];
        dx_local[i+1] += tmp[2];
        dy_local[i]   += tmp[1];
        dy_local[i+1] += tmp[3];
      } //end of big double loop
      
      #pragma omp single
      for( ; i < nodes; ++i ) {
        // clear displacement vectors
        dxi = 0.;
        dyi = 0.;
        
        pxi = position[i*dim];
        pyi = position[i*dim+1];
        
        for( j = i+1; j < nodes; ++j ) {
          dx = pxi - position[j*dim];
          dy = pyi - position[j*dim+1];
          len = dx*dx + dy*dy;
          
          // avoid division by zero
          if( len == 0. ) {
            dx = rand() * 1e-9;
            dy = rand() * 1e-9;
            len = dx*dx + dy*dy;
          }
          
          dxi += dx/len;
          dyi += dy/len;
          dx_local[j] -= dx/len;
          dy_local[j] -= dy/len;   
        }
        dx_local[i] += dxi;
        dy_local[i] += dyi;
      }
      
      #pragma omp for
      for(i = 0; i < nodes; ++i ) {

        dxi = dplx[i];
        dxi += dx_local[i];
        dyi = dply[i];
        dyi += dy_local[i];
        
        px  = position[i*dim];
        py  = position[i*dim+1];

        dxi += rand() * 1e-9;
        dyi += rand() * 1e-9;
        len = sqrt( dxi*dxi + dyi*dyi );
        
        real_dx = ( fabs(dxi) < temp ) ? dxi : temp;
        real_dy = ( fabs(dyi) < temp ) ? dyi : temp;
        
        // avoid division by zero
        if( len > 0. ) {
          px += (dxi / len) * real_dx;
          py += (dyi / len) * real_dy;
          
          position[i*dim] = px;
          position[i*dim+1] = py;
        }
      }
      aligned_free( tmp );
    } // end of parallel section
    temp -= difftemp;
  } // end of iteration loop
  
  // free allocated memory
  free( dplx );
  free( dply );
  free( space );
  
  gettimeofday(&end, NULL);
  double time_elapsed = gettime(start,end);
  printf("Layouting this graph took %1.4f sec\n", time_elapsed);
}

/* layout fruchterman reingold
* input: 
* - problem size: #nodes * dimension
* - edge list
*   2d array, first dimension source, second dim dest
* - MAX_IT:  number of iterations
* - temp:    some sort of temperature 
*/
// NOTE:
// - idea: use 1D array for easier vectoriztion afterwards.
// - should also reduce cachline trashing
// is about 5% faster than float baseline
void layout_fr_stride( float *position, int nodes, 
                       int *es, int edges, int dummy_dim, 
                       int MAX_IT, float temp ) {
  int n, i, j;
  int source, target;
  
  float dx, dy, len;
  float *dplx = (float*) aligned_malloc(nodes*sizeof(float));
  float *dply = (float*) aligned_malloc(nodes*sizeof(float));
  
  float real_dx, real_dy;
  float difftemp = temp / MAX_IT;
  
  // first, set random positions
  const unsigned int dim = 2;
  nodes /= dim;
  initialize( position, nodes, dim );
  
  float *xpos = position;
  float *ypos = position + nodes;
  
  // do this for some iterations
  for( n=0; n < MAX_IT; ++n ) {
    
    // clear displacement vectors
    for( i=0; i < nodes; ++i ) {
      dplx[i] = 0.;
      dply[i] = 0.;
    }
    
    // repulsive forces
    for( i = 0; i < nodes; ++i ) {
      for( j = i+1; j < nodes; ++j ) {
        dx = xpos[i] - xpos[j];
        dy = ypos[i] - ypos[j];
        len = dx*dx + dy*dy;
        
        // avoid division by zero
        if( len == 0. ) {
          dx = rand() * 1e-9;
          dy = rand() * 1e-9;
          len = dx*dx + dy*dy;
        }
        
        dplx[i] += dx/len;
        dply[i] += dy/len;
        dplx[j] -= dx/len;
        dply[j] -= dy/len;
      }
    }
    
    // attractive forces
    for( i = 0; i < edges; ++i ) {
      source = es[i*dim];
      target = es[i*dim+1];
      
      dx = xpos[source] - xpos[target];
      dy = ypos[source] - ypos[target];
      len = sqrt(dx*dx + dy*dy);
      
      dplx[source] -= dx * len;
      dply[source] -= dy * len;
      dplx[target] += dx * len;
      dply[target] += dy * len;
    }
    
    // finally update the positions
    for( i=0; i < nodes; ++i ) {
      dx = dplx[i] + rand() * 1e-9;
      dy = dply[i] + rand() * 1e-9;
      len = sqrt( dx*dx + dy*dy );
      
      real_dx = ( fabs(dx) < temp ) ? dx : temp;
      real_dy = ( fabs(dy) < temp ) ? dy : temp;
      
      // avoid division by zero
      if( len > 0. ) {
        xpos[i] += (dx / len) * real_dx;
        ypos[i] += (dy / len) * real_dy;
      }
    }
    
    temp -= difftemp;
  } // end of iteration loop
  
  // free allocated memory
  aligned_free( dplx );
  aligned_free( dply );
}

/* layout fruchterman reingold
* input: 
* - problem size: #nodes * dimension
* - edge list
*   2d array, first dimension source, second dim dest
* - MAX_IT:  number of iterations
* - temp:    some sort of temperature 
*/
// NOTE:
// version building on layout_fr_stride
// this version has a speedup over layout_fr_stride of 2.7 and is 
// therefore way better than anything else
void layout_fr_stride_simd( float *position, int nodes, 
                       int *es, int edges, int dummy_dim, 
                       int MAX_IT, float temp ) {
  int n, i, j;
  int source, target;
  
  float *dplx = (float*) aligned_malloc(nodes*sizeof(float));
  float *dply = (float*) aligned_malloc(nodes*sizeof(float));
  
  float *spacex = (float*) aligned_malloc(4*sizeof(float));
  float *spacey = (float*) aligned_malloc(4*sizeof(float));
  
  float dx, dy, len;
  float real_dx, real_dy;
  float difftemp = temp / MAX_IT;
  
  // first, set random positions
  const unsigned int dim = 2;
  nodes /= dim;
  initialize( position, nodes, dim );
  
  float *xpos = position;
  float *ypos = position + nodes;
  
  __m128 zero = _mm_setzero_ps();
  
  // do this for some iterations
  for( n=0; n < MAX_IT; ++n ) {
    
    // clear displacement vectors
    for( i=0; i < nodes; ++i ) {
      dplx[i] = 0.;
      dply[i] = 0.;
    }
    
    // repulsive forces
    for( i = 0; i < nodes-3; i+=4 ) {
      // reset displacement values for i
      __m128 dplxi = _mm_setzero_ps();
      __m128 dplyi = _mm_setzero_ps();

      // load i values
      __m128 xi = _mm_loadu_ps( &xpos[i] );      
      __m128 yi = _mm_loadu_ps( &ypos[i] );

      for( j = i+1; j < nodes; ++j ) {
        // load j values
        __m128 xj = _mm_load1_ps( &xpos[j] );
        __m128 yj = _mm_load1_ps( &ypos[j] );
        
        // define distances and length
        __m128 dx = _mm_sub_ps( xi, xj );
        __m128 dy = _mm_sub_ps( yi, yj );
        __m128 len = _mm_add_ps( _mm_mul_ps(dx, dx), _mm_mul_ps(dy, dy) );
        
        // NOTE
        // this is necessary to avoid division by zero
        // but it is very costly, so do it only if necessary
        if( _mm_movemask_ps( _mm_cmpeq_ps(zero, len) ) != 0 ) {
          _mm_store_ps( spacex, len );
          if( spacex[0] == 0. ) {
            float dx = rand() * 1e-9;
            float dy = rand() * 1e-9;
            spacex[0] = dx*dx + dy*dy;
          }
          if( spacex[1] == 0. ) {
            float dx = rand() * 1e-9;
            float dy = rand() * 1e-9;
            spacex[1] = dx*dx + dy*dy;
          }
          if( spacex[2] == 0. ) {
            float dx = rand() * 1e-9;
            float dy = rand() * 1e-9;
            spacex[2] = dx*dx + dy*dy;
          }
          if( spacex[3] == 0. ) {
            float dx = rand() * 1e-9;
            float dy = rand() * 1e-9;
            spacex[3] = dx*dx + dy*dy;
          }
          len = _mm_load_ps( spacex );
        }
       
        __m128 tmpx = _mm_div_ps( dx, len );
        __m128 tmpy = _mm_div_ps( dy, len );
        
        // update displacement vectors for i
        dplxi = _mm_add_ps( dplxi, tmpx );
        dplyi = _mm_add_ps( dplyi, tmpy );

        // update displacement vectors for j
        _mm_store_ps( spacex, tmpx );
        _mm_store_ps( spacey, tmpy );
        float sumx = spacex[0] + spacex[1];
        float sumx2 = spacex[2] + spacex[3];
        dplx[j] -= (sumx + sumx2);
        float sumy = spacey[0] + spacey[1];
        float sumy2 = spacey[2] + spacey[3];
        dply[j] -= (sumy + sumy2);
      }
      __m128 tmpx = _mm_load_ps( &dplx[i] );
      __m128 tmpy = _mm_load_ps( &dply[i] );
      
      dplxi = _mm_add_ps( tmpx, dplxi );
      dplyi = _mm_add_ps( tmpy, dplyi );
      
      _mm_store_ps( &dplx[i], dplxi );
      _mm_store_ps( &dply[i], dplyi );
    }
    
    // repulsive forces
    for( ; i < nodes; ++i ) {
      for( j = i+1; j < nodes; ++j ) {
        dx = xpos[i] - xpos[j];
        dy = ypos[i] - ypos[j];
        len = dx*dx + dy*dy;
        
        // avoid division by zero
        if( len == 0. ) {
          dx = rand() * 1e-9;
          dy = rand() * 1e-9;
          len = dx*dx + dy*dy;
        }
        
        dplx[i] += dx/len;
        dply[i] += dy/len;
        dplx[j] -= dx/len;
        dply[j] -= dy/len;
      }
    }
    
    // attractive forces
    for( i = 0; i < edges; ++i ) {
      source = es[i*dim];
      target = es[i*dim+1];
      
      dx = xpos[source] - xpos[target];
      dy = ypos[source] - ypos[target];
      len = sqrt(dx*dx + dy*dy);
      
      dplx[source] -= dx * len;
      dply[source] -= dy * len;
      dplx[target] += dx * len;
      dply[target] += dy * len;
    }
    
    // finally update the positions
    for( i=0; i < nodes; ++i ) {
      dx = dplx[i] + rand() * 1e-9;
      dy = dply[i] + rand() * 1e-9;
      len = sqrt( dx*dx + dy*dy );
      
      real_dx = ( fabs(dx) < temp ) ? dx : temp;
      real_dy = ( fabs(dy) < temp ) ? dy : temp;
      
      // avoid division by zero
      if( len > 0. ) {
        xpos[i] += (dx / len) * real_dx;
        ypos[i] += (dy / len) * real_dy;
      }
    }
    
    temp -= difftemp;
  } // end of iteration loop
  
  // free allocated memory
  aligned_free( dplx );
  aligned_free( dply );
  aligned_free( spacex );
  aligned_free( spacey );
}

/* layout fruchterman reingold
* input: 
* - problem size: #nodes * dimension
* - edge list
*   2d array, first dimension source, second dim dest
* - MAX_IT:  number of iterations
* - temp:    some sort of temperature 
*/
// NOTE:
// - idea: use 1D array for easier vectoriztion afterwards.
// - should also reduce cachline trashing
// is about 5% faster than float baseline
void layout_fr_stride_omp( float *position, int nodes, 
                       int *es, int edges, int dummy_dim, 
                       int MAX_IT, float temp ) {
  struct timeval start, end;
  gettimeofday(&start, NULL);
  int n, i, j;
  int source, target;
  
  float dx, dy, len;
  float px, py, dxi, dyi;
  float *dplx = (float*) aligned_malloc(nodes*sizeof(float));
  float *dply = (float*) aligned_malloc(nodes*sizeof(float));
  float *space = (float*) aligned_malloc(2*nodes*sizeof(float));
  
  float real_dx, real_dy;
  float difftemp = temp / MAX_IT;
  
  // first, set random positions
  const unsigned int dim = 2;
  nodes /= dim;
  initialize( position, nodes, dim );
  
  float *xpos = position;
  float *ypos = position + nodes;
  
  // do this for some iterations
  for( n=0; n < MAX_IT; ++n ) {
    
    // clear displacement vectors
    for( i=0; i < nodes; ++i ) {
      dplx[i] = 0.;
      dply[i] = 0.;
    }
    
    // attractive forces
    for( i = 0; i < edges; ++i ) {
      source = es[i*dim];
      target = es[i*dim+1];
      
      dx = xpos[source] - xpos[target];
      dy = ypos[source] - ypos[target];
      len = sqrt(dx*dx + dy*dy);
      
      dplx[source] -= dx * len;
      dply[source] -= dy * len;
      dplx[target] += dx * len;
      dply[target] += dy * len;
    }
    
    // repulsive forces
    //
    // NOTE
    // private displacement vectors are needed to avoid race conditions
    #pragma omp parallel \
    shared( xpos, ypos, dplx, dply, nodes, temp) \
    private(i, j, dx, dy, len, px, py, dxi, dyi, real_dx, real_dy) 
    {
      float *dx_local = space;
      float *dy_local = space + nodes;
      for(i = 0; i < nodes; ++i ) {
        dx_local[i] = 0.;
        dy_local[i] = 0.;
      }
      
      #pragma omp for schedule(dynamic)
      for( i = 0; i < nodes; ++i ) {
        dxi = 0.;
        dyi = 0.;
        
        px = xpos[i];
        py = ypos[i];
        for( j = i+1; j < nodes; ++j ) {
          dx = px - xpos[j];
          dy = py - ypos[j];
          len = dx*dx + dy*dy;
          
          // avoid division by zero
          if( len == 0. ) {
            dx = rand() * 1e-9;
            dy = rand() * 1e-9;
            len = dx*dx + dy*dy;
          }
          
          dxi += dx/len;
          dyi += dy/len;
          dx_local[j] -= dx/len;
          dy_local[j] -= dy/len;
        }
        dx_local[i] += dxi;
        dy_local[i] += dyi;
      } // end of double loop
      
      // finally update the positions
      #pragma omp for
      for( i=0; i < nodes; ++i ) {
        dx = dplx[i];
        dx += dx_local[i];
        dx += rand() * 1e-9;
        dy = dply[i];
        dy += dy_local[i];
        dy += rand() * 1e-9;
        
        len = sqrt( dx*dx + dy*dy );
        
        real_dx = ( fabs(dx) < temp ) ? dx : temp;
        real_dy = ( fabs(dy) < temp ) ? dy : temp;
        
        // avoid division by zero
        if( len > 0. ) {
          xpos[i] += (dx / len) * real_dx;
          ypos[i] += (dy / len) * real_dy;
        }
      }
    } // end of parallel section
    
    temp -= difftemp;
  } // end of iteration loop
  
  // free allocated memory
  aligned_free( dplx );
  aligned_free( dply );
  aligned_free( space );
  
  gettimeofday(&end, NULL);
  double time_elapsed = gettime(start,end);
  printf("Layouting this graph took %1.4f sec\n", time_elapsed);
}

/* layout fruchterman reingold
* input: 
* - problem size: #nodes * dimension
* - edge list
*   2d array, first dimension source, second dim dest
* - MAX_IT:  number of iterations
* - temp:    some sort of temperature 
*/
// NOTE:
// - idea: use 1D array for easier vectoriztion afterwards.
// - should also reduce cachline trashing
// is about 5% faster than float baseline
void layout_fr_stride_omp_simd( float *position, int nodes, 
                       int *es, int edges, int dummy_dim, 
                       int MAX_IT, float temp ) {
  struct timeval start, end;
  gettimeofday(&start, NULL);
  int n, i, j;
  int source, target;
  
  float dx, dy, len;
  float px, py;
  float *dplx = (float*) aligned_malloc(nodes*sizeof(float));
  float *dply = (float*) aligned_malloc(nodes*sizeof(float));
  float *space = (float*) aligned_malloc(2*nodes*sizeof(float));
  
  float real_dx, real_dy;
  float difftemp = temp / MAX_IT;
  
  // first, set random positions
  const unsigned int dim = 2;
  nodes /= dim;
  initialize( position, nodes, dim );
  
  float *xpos = position;
  float *ypos = position + nodes;
  
  __m128 zero = _mm_setzero_ps();
  
  // do this for some iterations
  for( n=0; n < MAX_IT; ++n ) {
    
    // clear displacement vectors
    for( i=0; i < nodes; ++i ) {
      dplx[i] = 0.;
      dply[i] = 0.;
    }
    
    // attractive forces
    for( i = 0; i < edges; ++i ) {
      source = es[i*dim];
      target = es[i*dim+1];
      
      dx = xpos[source] - xpos[target];
      dy = ypos[source] - ypos[target];
      len = sqrt(dx*dx + dy*dy);
      
      dplx[source] -= dx * len;
      dply[source] -= dy * len;
      dplx[target] += dx * len;
      dply[target] += dy * len;
    }
    
    // repulsive forces
    //
    // NOTE
    // private displacement vectors are needed to avoid race conditions
    #pragma omp parallel \
    shared( xpos, ypos, dplx, dply, nodes, temp) \
    private(i, j, dx, dy, len, px, py, real_dx, real_dy) 
    {
      float *spacex = (float*) aligned_malloc(4*sizeof(float));
      float *spacey = (float*) aligned_malloc(4*sizeof(float));
      
      #pragma omp for schedule(dynamic)
      for( i = 0; i < nodes-3; i+=4 ) {
        // reset displacement values for i
        __m128 dplxi = _mm_setzero_ps();
        __m128 dplyi = _mm_setzero_ps();

        // load i values
        __m128 xi = _mm_loadu_ps( &xpos[i] );      
        __m128 yi = _mm_loadu_ps( &ypos[i] );

        for( j = i+1; j < nodes; ++j ) {
          // load j values
          __m128 xj = _mm_load1_ps( &xpos[j] );
          __m128 yj = _mm_load1_ps( &ypos[j] );
          
          // define distances and length
          __m128 dx = _mm_sub_ps( xi, xj );
          __m128 dy = _mm_sub_ps( yi, yj );
          __m128 len = _mm_add_ps( _mm_mul_ps(dx, dx), _mm_mul_ps(dy, dy) );
          
          // NOTE
          // this is necessary to avoid division by zero
          // but it is very costly, so do it only if necessary
          if( _mm_movemask_ps( _mm_cmpeq_ps(zero, len) ) != 0 ) {
            _mm_store_ps( spacex, len );
            if( spacex[0] == 0. ) {
              float dx = rand() * 1e-9;
              float dy = rand() * 1e-9;
              spacex[0] = dx*dx + dy*dy;
            }
            if( spacex[1] == 0. ) {
              float dx = rand() * 1e-9;
              float dy = rand() * 1e-9;
              spacex[1] = dx*dx + dy*dy;
            }
            if( spacex[2] == 0. ) {
              float dx = rand() * 1e-9;
              float dy = rand() * 1e-9;
              spacex[2] = dx*dx + dy*dy;
            }
            if( spacex[3] == 0. ) {
              float dx = rand() * 1e-9;
              float dy = rand() * 1e-9;
              spacex[3] = dx*dx + dy*dy;
            }
            len = _mm_load_ps( spacex );
          }
        
          __m128 tmpx = _mm_div_ps( dx, len );
          __m128 tmpy = _mm_div_ps( dy, len );
          
          // update displacement vectors for i
          dplxi = _mm_add_ps( dplxi, tmpx );
          dplyi = _mm_add_ps( dplyi, tmpy );

          // update displacement vectors for j
          _mm_store_ps( spacex, tmpx );
          _mm_store_ps( spacey, tmpy );
          float sumx = spacex[0] + spacex[1];
          float sumx2 = spacex[2] + spacex[3];
          dplx[j] -= (sumx + sumx2);
          float sumy = spacey[0] + spacey[1];
          float sumy2 = spacey[2] + spacey[3];
          dply[j] -= (sumy + sumy2);
        }
        __m128 tmpx = _mm_load_ps( &dplx[i] );
        __m128 tmpy = _mm_load_ps( &dply[i] );
        
        dplxi = _mm_add_ps( tmpx, dplxi );
        dplyi = _mm_add_ps( tmpy, dplyi );
        
        _mm_store_ps( &dplx[i], dplxi );
        _mm_store_ps( &dply[i], dplyi );
      }
      
      // repulsive forces
      #pragma omp single
      for(i = nodes/4*4 ; i < nodes; ++i ) {
        px = xpos[i];
        py = ypos[i];
        for( j = i+1; j < nodes; ++j ) {
          dx = px - xpos[j];
          dy = py - ypos[j];
          len = dx*dx + dy*dy;
          
          // avoid division by zero
          if( len == 0. ) {
            dx = rand() * 1e-9;
            dy = rand() * 1e-9;
            len = dx*dx + dy*dy;
          }
          
          dplx[i] += dx/len;
          dply[i] += dy/len;
          dplx[j] -= dx/len;
          dply[j] -= dy/len;
        }
      }
      
      // finally update the positions
      #pragma omp for
      for( i=0; i < nodes; ++i ) {
        dx = dplx[i] + rand() * 1e-9;
        dy = dply[i] + rand() * 1e-9;
        
        len = sqrt( dx*dx + dy*dy );
        
        real_dx = ( fabs(dx) < temp ) ? dx : temp;
        real_dy = ( fabs(dy) < temp ) ? dy : temp;
        
        // avoid division by zero
        if( len > 0. ) {
          xpos[i] += (dx / len) * real_dx;
          ypos[i] += (dy / len) * real_dy;
        }
      }
      
      aligned_free( spacex );
      aligned_free( spacey );
    } // end of parallel section
    
    temp -= difftemp;
  } // end of iteration loop
  
  // free allocated memory
  aligned_free( dplx );
  aligned_free( dply );
  aligned_free( space );
  
  gettimeofday(&end, NULL);
  double time_elapsed = gettime(start,end);
  printf("Layouting this graph took %1.4f sec\n", time_elapsed);
}
