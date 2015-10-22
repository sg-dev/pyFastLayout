#include "fastlayout.h"
// aligned memory
#include "memory.c"

inline double gettime( struct timeval start, struct timeval end ) {
  return (double)((end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)/1e6);
}

/* layout fruchterman reingold
* input: 
* - problem size: #nodes * dimension
* - edge list
*   2d array, first dimension source, second dim dest
* - MAX_IT:  number of iterations
* - temp:    some sort of temperature 
*/
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
  float sqrt_nodes = sqrt(nodes);
  for( i = 0; i < nodes*dim; ++i ) {
    // math.sqrt(nodes) * np.random.random_sample( (nodes, dimension) ) - math.sqrt(nodes)/2
    position[i] = sqrt_nodes * (rand() / (float) ((unsigned)RAND_MAX+1)) - sqrt_nodes/2;
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
  float sqrt_nodes = sqrt(nodes);
  for( i = 0; i < nodes*dim; ++i ) {
    // math.sqrt(nodes) * np.random.random_sample( (nodes, dimension) ) - math.sqrt(nodes)/2
    position[i] = sqrt_nodes * (rand() / (float) ((unsigned)RAND_MAX+1)) - sqrt_nodes/2;
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
// NOTE
// - version preparing 4 way float vetorization
// - idea: unrolling iloop by 2 and jloop by 2
void layout_fr_2x2( float *position, int nodes, 
                    int *es, int edges, int dummy_dim, 
                    int MAX_IT, float temp ) {
  int n, i, j;
  int source, target;
  
  float dx, dy, len;
  float dx2, dy2, len2;
  float dx3, dy3, len3, dx4, dy4, len4;
  float *dplx = (float*) aligned_malloc(nodes*sizeof(float));
  float *dply = (float*) aligned_malloc(nodes*sizeof(float));
  
  float real_dx, real_dy;
  float difftemp = temp / MAX_IT;
  float pxi, pyi, dplxi, dplyi;
  
  // first, set random positions
  const unsigned int dim = 2;
  nodes /= dim;
  float sqrt_nodes = sqrt(nodes);
  for( i = 0; i < nodes*dim; ++i ) {
    // math.sqrt(nodes) * np.random.random_sample( (nodes, dimension) ) - math.sqrt(nodes)/2
    position[i] = sqrt_nodes * (rand() / (float) ((unsigned)RAND_MAX+1)) - sqrt_nodes/2;
  }
  
  // do this for some iterations
  for( n=0; n < MAX_IT; ++n ) {
    // clear displacement vectors
    for( i = 0; i < nodes; ++i ) {
      dplx[i] = 0;
      dply[i] = 0;
    }
    
    // repulsive forces
    for( i = 0; i < nodes; ++i ) {
      pxi = position[i*dim];
      pyi = position[i*dim+1];
      
      dplxi = 0;
      dplyi = 0;
      
      // do all j's here, exept possibly the last one
      for( j = i+1; (j+1) < nodes; j+=4 ) {
        dx = pxi - position[j*dim];
        dy = pyi - position[j*dim+1];
        
        dx2 = pxi - position[(j+1)*dim];
        dy2 = pyi - position[(j+1)*dim+1];
        
        dx3 = pxi - position[(j+2)*dim];
        dy3 = pyi - position[(j+2)*dim+1];
        
        dx4 = pxi - position[(j+3)*dim];
        dy4 = pyi - position[(j+3)*dim+1];
        
        len = dx*dx + dy*dy;
        len2 = dx2*dx2 + dy2*dy2;
        len3 = dx3*dx3 + dy3*dy3;
        len4 = dx4*dx4 + dy4*dy4;
        
        // avoid division by zero
        if( len == 0. ) {
          dx = rand() * 1e-9;
          dy = rand() * 1e-9;
          len = dx*dx + dy*dy;
        }
        if( len2 == 0. ) {
          dx = rand() * 1e-9;
          dy = rand() * 1e-9;
          len2 = dx*dx + dy*dy;
        }
        if( len3 == 0. ) {
          dx = rand() * 1e-9;
          dy = rand() * 1e-9;
          len3 = dx*dx + dy*dy;
        }
        if( len4 == 0. ) {
          dx = rand() * 1e-9;
          dy = rand() * 1e-9;
          len4 = dx*dx + dy*dy;
        }
        
        dplxi += (dx/len + dx2/len2 + dx3/len3 + dx4/len4);
        dplyi += (dy/len + dy2/len2 + dy3/len3 + dy4/len4);
        dplx[j]   -= dx/len;
        dplx[j+1] -= dx2/len2;
        dplx[j+2] -= dx3/len3;
        dplx[j+3] -= dx4/len4;
        dply[j]   -= dy/len;
        dply[j+1] -= dy2/len2;
        dply[j+2] -= dy3/len3;
        dply[j+3] -= dy4/len4;
      }
      // and the last one is again separte
      for( ; j < nodes; ++j ) {
        dx = pxi - position[j*dim];
        dy = pyi - position[j*dim+1];
        len = dx*dx + dy*dy;
        
        // avoid division by zero
        if( len == 0. ) {
          dx = rand() * 1e-9;
          dy = rand() * 1e-9;
          len = dx*dx + dy*dy;
        }
        
        dplxi += dx/len;
        dplyi += dy/len;
        dplx[j] -= dx/len;
        dply[j] -= dy/len;
      }
      
      dplx[i] += dplxi;
      dply[i] += dplyi;
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
  float sqrt_nodes = sqrt(nodes);
  for( i = 0; i < nodes*dim; ++i ) {
    // math.sqrt(nodes) * np.random.random_sample( (nodes, dimension) ) - math.sqrt(nodes)/2
    position[i] = sqrt_nodes * (rand() / (float) ((unsigned)RAND_MAX+1)) - sqrt_nodes/2;
  }
  
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
void layout_fr_loop_reorder( float *position, int nodes,
                             int *es, int edges, int dummy_dim,  
                             int MAX_IT, float temp ) {
  int n, i, j;
  int source, target;
  
  float dx, dy, len, dxi, dyi, px, py;
  float *dplx = (float*) malloc(nodes*sizeof(float));
  float *dply = (float*) malloc(nodes*sizeof(float));
  
  float real_dx, real_dy;
  float difftemp = temp / MAX_IT;
  
  // first, set random positions
  const unsigned int dim = 2;
  nodes /= dim;
  float sqrt_nodes = sqrt(nodes);
  for( i = 0; i < nodes*dim; ++i ) {
    // math.sqrt(nodes) * np.random.random_sample( (nodes, dimension) ) - math.sqrt(nodes)/2
    position[i] = sqrt_nodes * (rand() / (float) ((unsigned)RAND_MAX+1)) - sqrt_nodes/2;
  }
  
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
    for( i = 0; i < nodes; ++i ) {
      // clear displacement vectors
      dxi = dplx[i];
      dyi = dply[i];
      
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
        dplx[j] -= dx/len;
        dply[j] -= dy/len;
      }

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
    
    temp -= difftemp;
  } // end of iteration loop
  
  // free allocated memory
  free( dplx );
  free( dply );
}

/* layout fruchterman reingold
* input: 
* - problem size: #nodes * dimension
* - edge list
*   2d array, first dimension source, second dim dest
* - MAX_IT:  number of iterations
* - temp:    some sort of temperature 
*/
void layout_fr_block_iloop2( float *position, int nodes, 
                             int *es, int edges, int dummy_dim, 
                             int MAX_IT, float temp ){
  int n, i, j;
  int source, target;
  
  float dx, dy, len, dplxi, dplyi, px, py;
  float dxii, dyii, lenii, dplxii, dplyii, pxii, pyii;
  float pxj, pyj;
  float *dplx = (float*) aligned_malloc(nodes*sizeof(float));
  float *dply = (float*) aligned_malloc(nodes*sizeof(float));
  
  float real_dx, real_dy;
  float difftemp = temp / MAX_IT;
  
  // first, set random positions
  const unsigned int dim = 2;
  nodes /= dim;
  float sqrt_nodes = sqrt(nodes);
  for( i = 0; i < nodes*dim; ++i ) {
    // math.sqrt(nodes) * np.random.random_sample( (nodes, dimension) ) - math.sqrt(nodes)/2
    position[i] = sqrt_nodes * (rand() / (float) ((unsigned)RAND_MAX+1)) - sqrt_nodes/2;
  }
  
  // do this for some iterations
  for( n=0; n < MAX_IT; ++n ) {
    
    // repulsive forces
    for( i = 0; (i+1) < nodes; i+=2 ) {
      // clear displacement vectors
      dplxi = 0.;
      dplyi = 0.;
      dplxii = 0.;
      dplyii = 0.;
      
      px = position[i*dim];
      py = position[i*dim+1];
      pxii = position[(i+1)*dim];
      pyii = position[(i+1)*dim+1];
      
      // NOTE:
      // j=i+1 needs special care
      // p_j needs not to be process with p_{i+1} in this case!
      // see also simd versions
      j = i+1;
      if( j < nodes ) {
        pxj = position[j*dim];
        pyj = position[j*dim+1];
        
        dx = px - pxj;
        dy = py - pyj;
        len = dx*dx + dy*dy;
        
        // avoid division by zero
        if( len == 0. ) {
          dx = rand() * 1e-9;
          dy = rand() * 1e-9;
          len = dx*dx + dy*dy;
        }
        
        dplxi += dx/len;
        dplyi += dy/len;
        dplx[j] -= dx/len;
        dply[j] -= dy/len; 
      }
      // do all other j's now
      for( j = i+2; j < nodes; ++j ) {
        pxj = position[j*dim];
        pyj = position[j*dim+1];
        
        dx = px - pxj;
        dy = py - pyj;
        dxii = pxii - pxj;
        dyii = pyii - pyj;
        len = dx*dx + dy*dy;
        lenii = dxii*dxii + dyii*dyii;
        
        // avoid division by zero
        if( len == 0. ) {
          dx = rand() * 1e-9;
          dy = rand() * 1e-9;
          len = dx*dx + dy*dy;
        }
        if(lenii == 0.) {
          dxii = rand() * 1e-9;
          dyii = rand() * 1e-9;
          lenii = dxii*dxii + dyii*dyii;
        }
        
        dplxi += dx/len;
        dplyi += dy/len;
        dplxii += dxii/lenii;
        dplyii += dyii/lenii;
        dplx[j] -= (dx/len + dxii/lenii);
        dply[j] -= (dy/len + dyii/lenii); 
      }
      dplx[i] = dplxi;
      dply[i] = dplyi;
      dplx[i+1] = dplxii;
      dply[i+1] = dplyii;
    }
    
    for( ; i < nodes; ++i ) {
      // clear displacement vectors
      dplxi = 0.;
      dplyi = 0.;
      
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
        
        dplxi += dx/len;
        dplyi += dy/len;
        dplx[j] -= dx/len;
        dply[j] -= dy/len;   
      }
      dplx[i] = dplxi;
      dply[i] = dplyi;
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



// ------------------------------------------------------------------------
// OPENMP versions
// ------------------------------------------------------------------------  

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
  const int dim = 2;
  nodes /= dim;
  const float sqrt_nodes = sqrt(nodes);
  for( i = 0; i < nodes*dim; ++i ) {
    // math.sqrt(nodes) * np.random.random_sample( (nodes, dimension) ) - math.sqrt(nodes)/2
    position[i] = sqrt_nodes * (rand() / (float) ((unsigned)RAND_MAX+1)) - sqrt_nodes/2;
  }
  
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
  const int dim = 2;
  nodes /= dim;
  const float sqrt_nodes = sqrt(nodes);
  for( i = 0; i < nodes*dim; ++i ) {
    // math.sqrt(nodes) * np.random.random_sample( (nodes, dimension) ) - math.sqrt(nodes)/2
    position[i] = sqrt_nodes * (rand() / (float) ((unsigned)RAND_MAX+1)) - sqrt_nodes/2;
  }
  
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



//-------------------------------------------------------------------------
// SEE versions
//-------------------------------------------------------------------------
// #ifdef SEE

/* layout fruchterman reingold
* 
* assumption:
* 1) area to layout is [#nodes] x [#nodes]
* 2) positions hold positions that are uniformly random in [-#nodes/2:+#nodes/2]
* 
* input: 
* - random positions:
*   2d array, length: no_nodes, width: dimension = 2 
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
  float sqrt_nodes = sqrt(nodes);
  for( i = 0; i < nodes*dim; ++i ) {
    // math.sqrt(nodes) * np.random.random_sample( (nodes, dimension) ) - math.sqrt(nodes)/2
    position[i] = sqrt_nodes * (rand() / (float) ((unsigned)RAND_MAX+1)) - sqrt_nodes/2;
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

/* layout fruchterman reingold
* 
* assumption:
* 1) area to layout is [#nodes] x [#nodes]
* 2) positions hold positions that are uniformly random in [-#nodes/2:+#nodes/2]
* 
* input: 
* - random positions:
*   2d array, length: no_nodes, width: dimension = 2 
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
  float sqrt_nodes = sqrt(nodes);
  for( i = 0; i < nodes*dim; ++i ) {
    // math.sqrt(nodes) * np.random.random_sample( (nodes, dimension) ) - math.sqrt(nodes)/2
    position[i] = sqrt_nodes * (rand() / (float) ((unsigned)RAND_MAX+1)) - sqrt_nodes/2;
  }
  
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
* 
* assumption:
* 1) area to layout is [#nodes] x [#nodes]
* 2) positions hold positions that are uniformly random in [-#nodes/2:+#nodes/2]
* 
* input: 
* - random positions:
*   2d array, length: no_nodes, width: dimension = 2 
* - edge list
*   2d array, first dimension source, second dim dest
* - MAX_IT:  number of iterations
* - temp:    some sort of temperature 
*/
// TODO
// not yet working 4way float sse version
// and also not faster than layout_fr_sse_float
void layout_fr_simd_jloop( float *position, int nodes,
                           int *es, int edges, int dummy_dim, 
                           int MAX_IT, float temp ) {
  int n, i, j;
  int source, target;
  
  float dx, dy, len, len2;
  float dxi, dyi;
  float *dplx = (float*) aligned_malloc(nodes*sizeof(float));
  float *dply = (float*) aligned_malloc(nodes*sizeof(float));
  float *tmp = (float*) aligned_malloc(4*sizeof(float));
  for( i = 0; i < nodes; ++i ) {
    dplx[i] = 0.;
    dply[i] = 0.;
  }
  
  float real_dx, real_dy;
  float difftemp = temp / MAX_IT;
  
  __m128 di, pi, pj, d, l;
  pi = _mm_setzero_ps();
  
  // first, set random positions
  const unsigned int dim = 2;
  nodes /= dim;
  float sqrt_nodes = sqrt(nodes);
  for( i = 0; i < nodes*dim; ++i ) {
    // math.sqrt(nodes) * np.random.random_sample( (nodes, dimension) ) - math.sqrt(nodes)/2
    position[i] = sqrt_nodes * (rand() / (float) ((unsigned)RAND_MAX+1)) - sqrt_nodes/2;
  }
  
  // do this for some iterations
  for( n=0; n < MAX_IT; ++n ) {
    
    for( i = 0; i < nodes; ++i ) {
      dplx[i] = 0;
      dply[i] = 0;
    }
    
    // repulsive forces
    for( i = 0; i < nodes; ++i ) {
      // clear displacement vectors
      di = _mm_setzero_ps();
      // load low  pi position[i*dim]
      // load high pi position[i*dim]
      _mm_loadl_pi(pi, (__m64*)&position[i*dim]);
      _mm_loadh_pi(pi, (__m64*)&position[i*dim]);
      
      for( j = i+1; (j+1) < nodes; j+=2 ) {
        pj = _mm_loadu_ps( &position[j*dim] );
        d = _mm_sub_ps( pi, pj );
        
        // store d*d in tmp
        _mm_store_ps(tmp, _mm_mul_ps(d, d));
        len = tmp[0] + tmp[1];
        len2 = tmp[2] + tmp[3];
        
        // avoid division by zero
        if( len == 0. ) {
          dx = rand() * 1e-9;
          dy = rand() * 1e-9;
          len = dx*dx + dy*dy;
        }
        if( len2 == 0. ) {
          dx = rand() * 1e-9;
          dy = rand() * 1e-9;
          len2 = dx*dx + dy*dy;
        }
        
        //          LSB
        //          |
        // load l = [len, len, len2, len2]
        l = _mm_set_ps(len2, len2, len, len);
        // d = d / l
        d = _mm_div_ps( d, l );
        // di = di + d
        di = _mm_add_ps(di, d);
        
        // store d to tmp
        _mm_store_ps(tmp, d);
        dplx[j]   -= tmp[0];
        dply[j]   -= tmp[1];
        dplx[j+1] -= tmp[2];
        dply[j+1] -= tmp[3];
      }
      // store di to tmp
      _mm_store_ps(tmp, di);
      dxi = tmp[0] + tmp[2];
      dyi = tmp[1] + tmp[3];
      
      for( ; j < nodes; ++j ) {
        dx = position[i*dim] - position[j*dim];
        dy = position[i*dim+1] - position[j*dim+1];
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
* 
* assumption:
* 1) area to layout is [#nodes] x [#nodes]
* 2) positions hold positions that are uniformly random in [-#nodes/2:+#nodes/2]
* 
* input: 
* - random positions:
*   2d array, length: no_nodes, width: dimension = 2 
* - edge list
*   2d array, first dimension source, second dim dest
* - MAX_IT:  number of iterations
* - temp:    some sort of temperature 
*/
// NOTE
// - version for 4 way float vetorization
// - idea: jloop by 4
void layout_fr_simd_1x4( float *position, int nodes,
                         int *es, int edges, int dummy_dim, 
                         int MAX_IT, float temp ) {
  int n, i, j;
  int source, target;
  
  float dx, dy, len;
  float *dplx = (float*) aligned_malloc(nodes*sizeof(float));
  float *dply = (float*) aligned_malloc(nodes*sizeof(float));
  float *tmp1 = (float*) aligned_malloc(4*sizeof(float));
  float *tmp2 = (float*) aligned_malloc(4*sizeof(float));
  
  float real_dx, real_dy;
  float difftemp = temp / MAX_IT;
  
  float pxi, pyi, dplxi, dplyi;
  
  __m128 pi, pj, pj2;
  __m128 d1, d2, lenv, lenv1, lenv2;
  __m128 t1, t2, dxv, dyv;
  
  // first, set random positions
  const unsigned int dim = 2;
  nodes /= dim;
  float sqrt_nodes = sqrt(nodes);
  for( i = 0; i < nodes*dim; ++i ) {
    // math.sqrt(nodes) * np.random.random_sample( (nodes, dimension) ) - math.sqrt(nodes)/2
    position[i] = sqrt_nodes * (rand() / (float) ((unsigned)RAND_MAX+1)) - sqrt_nodes/2;
  }
  
  // do this for some iterations
  for( n=0; n < MAX_IT; ++n ) {
    // clear displacement vectors
    for( i = 0; i < nodes; ++i ) {
      dplx[i] = 0;
      dply[i] = 0;
    }
    
    // repulsive forces
    for( i = 0; i < nodes; ++i ) {
      pxi = position[i*dim];
      pyi = position[i*dim+1];
      
      pi = _mm_loadl_pi( pi, (__m64*)&position[i*dim] );
      pi = _mm_loadh_pi( pi, (__m64*)&position[i*dim] );
      
      dplxi = 0;
      dplyi = 0;
      
      // do all j's here, exept possibly the last one
      for( j = i+1; (j+1) < nodes; j+=4 ) {
        pj = _mm_loadu_ps( &position[j*dim] );
        pj2 = _mm_loadu_ps( &position[(j+2)*dim] );
        
        dxv = _mm_loadu_ps( &dplx[j] );
        dyv = _mm_loadu_ps( &dply[j] );
        
        d1 = _mm_sub_ps( pi, pj );
        d2 = _mm_sub_ps( pi, pj2 );
        
        //       LSB
        //       |
        // len = [len, len2, len3, len4]
        lenv = _mm_hadd_ps( _mm_mul_ps(d1, d1), _mm_mul_ps(d2, d2) );
        
        // TODO:
        // avoid division by zero
        /*
        if( len == 0. ) {
          dx = rand() * 1e-9;
          dy = rand() * 1e-9;
          len = dx*dx + dy*dy;
        }
        if( len2 == 0. ) {
          dx = rand() * 1e-9;
          dy = rand() * 1e-9;
          len2 = dx*dx + dy*dy;
        }
        if( len3 == 0. ) {
          dx = rand() * 1e-9;
          dy = rand() * 1e-9;
          len3 = dx*dx + dy*dy;
        }
        if( len4 == 0. ) {
          dx = rand() * 1e-9;
          dy = rand() * 1e-9;
          len4 = dx*dx + dy*dy;
        }
        */
        
        // lenv1 = [len, len, len2, len2]
        // lenv2 = [len3, len3, len4, len4]
        lenv1 = _mm_shuffle_ps( lenv, lenv, _MM_SHUFFLE(1,1,0,0) );
        lenv2 = _mm_shuffle_ps( lenv, lenv, _MM_SHUFFLE(3,3,2,2) );
        
        // d1 = [  dx/len,   dy/len,  dx2/len2, dy2/len2 ]
        // d2 = [ dx3/len3, dy3/len3, dx4/len4, dy4/len4 ]
        d1 = _mm_div_ps( d1, lenv1 );
        d2 = _mm_div_ps( d2, lenv2 );
        t1 = _mm_shuffle_ps(d1, d2, _MM_SHUFFLE(2,0,2,0));
        t2 = _mm_shuffle_ps(d1, d2, _MM_SHUFFLE(3,1,3,1));
        dxv = _mm_sub_ps(dxv, t1);
        dyv = _mm_sub_ps(dyv, t2);
        
        // tmp1 = dx/len, dx2/len2, dx3/len3, dx4/len4
        _mm_store_ps( tmp1, t1 );
        _mm_store_ps( tmp2, t2 );
        _mm_storeu_ps( &dplx[j], dxv );
        _mm_storeu_ps( &dply[j], dyv );
        
        // store tmp1 and tmp2 at the right places
        dplxi += (tmp1[0] + tmp1[1] + tmp1[2] + tmp1[3]);
        dplyi += (tmp2[0] + tmp2[1] + tmp2[2] + tmp2[3]);
      }
      // the last one is separte
      for( ; j < nodes; ++j ) {
        dx = pxi - position[j*dim];
        dy = pyi - position[j*dim+1];
        len = dx*dx + dy*dy;
        
        // avoid division by zero
        if( len == 0. ) {
          dx = rand() * 1e-9;
          dy = rand() * 1e-9;
          len = dx*dx + dy*dy;
        }
        
        dplxi += dx/len;
        dplyi += dy/len;
        dplx[j] -= dx/len;
        dply[j] -= dy/len;
      }
      
      dplx[i] += dplxi;
      dply[i] += dplyi;
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
  aligned_free( tmp1 );
  aligned_free( tmp2 );
}

/* layout fruchterman reingold
* 
* assumption:
* 1) area to layout is [#nodes] x [#nodes]
* 2) positions hold positions that are uniformly random in [-#nodes/2:+#nodes/2]
* 
* input: 
* - random positions:
*   2d array, length: no_nodes, width: dimension = 2 
* - edge list
*   2d array, first dimension source, second dim dest
* - MAX_IT:  number of iterations
* - temp:    some sort of temperature 
*/
// TODO
// - not working but probably faster variant of the above function

/*
void layout_fr_float_1x4_sse( int *es, int edges, int dummy_dim, float *position, int nodes, int dim, int MAX_IT, float temp ) {
  struct timeval start, end;
  gettimeofday(&start, NULL);
  
  int n, i, j;
  int source, target;
  
  float dx, dy, len;
  float dx2, dy2, len2;
  float dx3, dy3, len3, dx4, dy4, len4;
  float *dplx = (float*) aligned_malloc(nodes*sizeof(float));
  float *dply = (float*) aligned_malloc(nodes*sizeof(float));
  float *tmp1 = (float*) aligned_malloc(4*sizeof(float));
  float *tmp2 = (float*) aligned_malloc(4*sizeof(float));
  
  float real_dx, real_dy;
  float difftemp = temp / MAX_IT;
  
  float pxi, pyi, dplxi, dplyi;
  
  __m128 pi, pj, pj2;
  __m128 d1, d2, lenv, lenv1, lenv2;
  __m128 savex, savey, dplxj, dplyj;
  
  // do this for some iterations
  for( n=0; n < MAX_IT; ++n ) {
    
    for( i = 0; i < n; ++i ) {
      dplx[i] = 0;
      dply[i] = 0;
    }
    
    // repulsive forces
    for( i = 0; i < nodes; ++i ) {
      pxi = position[i*dim];
      pyi = position[i*dim+1];
      
      pi = _mm_loadl_pi( pi, (__m64*)&position[i*dim] );
      pi = _mm_loadh_pi( pi, (__m64*)&position[i*dim] );
      
//       dplxi_v = _mm_setzero_ps();
//       dplyi_v = _mm_setzero_ps();
      
      dplxi = 0;
      dplyi = 0;
      
      // do all j's here, exept possibly the last one
      for( j = i+1; (j+1) < nodes; j+=4 ) {
        pj = _mm_loadu_ps( &position[j*dim] );
        pj2 = _mm_loadu_ps( &position[(j+2)*dim] );
        dplxj = _mm_loadu_ps( &dplx[j] );
        dplyj = _mm_loadu_ps( &dply[j] );
        
        d1 = _mm_sub_ps( pi, pj );
        d2 = _mm_sub_ps( pi, pj2 );
        
//         dx = pxi - position[j*dim];
//         dy = pyi - position[j*dim+1];
        
//         dx2 = pxi - position[(j+1)*dim];
//         dy2 = pyi - position[(j+1)*dim+1];
        
//         dx3 = pxi - position[(j+2)*dim];
//         dy3 = pyi - position[(j+2)*dim+1];
        
//         dx4 = pxi - position[(j+3)*dim];
//         dy4 = pyi - position[(j+3)*dim+1];
        
//         len = dx*dx + dy*dy;
//         len2 = dx2*dx2 + dy2*dy2;
//         len3 = dx3*dx3 + dy3*dy3;
//         len4 = dx4*dx4 + dy4*dy4;
        
        //       LSB
        //       |
        // len = [len, len2, len3, len4]
        lenv = _mm_hadd_ps( _mm_mul_ps(d1, d1), _mm_mul_ps(d2, d2) );
        
        // TODO:
        // avoid division by zero
//         if( len == 0. ) {
//           dx = rand() * 1e-9;
//           dy = rand() * 1e-9;
//           len = dx*dx + dy*dy;
//         }
//         if( len2 == 0. ) {
//           dx = rand() * 1e-9;
//           dy = rand() * 1e-9;
//           len2 = dx*dx + dy*dy;
//         }
//         if( len3 == 0. ) {
//           dx = rand() * 1e-9;
//           dy = rand() * 1e-9;
//           len3 = dx*dx + dy*dy;
//         }
//         if( len4 == 0. ) {
//           dx = rand() * 1e-9;
//           dy = rand() * 1e-9;
//           len4 = dx*dx + dy*dy;
//         }
        
        
        // lenv1 = [len, len, len2, len2]
        // lenv2 = [len3, len3, len4, len4]
        lenv1 = _mm_shuffle_ps( lenv, lenv, _MM_SHUFFLE(1,1,0,0) );
        lenv2 = _mm_shuffle_ps( lenv, lenv, _MM_SHUFFLE(3,3,2,2) );
        
        // d1 = [  dx/len,   dy/len,  dx2/len2, dy2/len2 ]
        // d2 = [ dx3/len3, dy3/len3, dx4/len4, dy4/len4 ]
        d1 = _mm_div_ps( d1, lenv1 );
        d2 = _mm_div_ps( d2, lenv2 );
        
        // savex = [  dx/len, dx2/len2, dx3/len3, dx4/len4 ]
        // savey = [  dy/len, dy2/len2, dy3/len3, dy3/len4 ]
        savex = _mm_shuffle_ps( d1, d2, _MM_SHUFFLE( 2, 0, 2, 0 ));
        savey = _mm_shuffle_ps( d1, d2, _MM_SHUFFLE( 3, 1, 3, 1 ));
        dplxj = _mm_sub_ps( dplxj, savex );
        dplyj = _mm_sub_ps( dplyj, savey );
        
        _mm_storeu_ps( &dplx[j], dplxj );
        _mm_storeu_ps( &dply[j], dplyj );
        _mm_store_ps( tmp1, savex );
        _mm_store_ps( tmp2, savey );
        

//         dplxi += (dx/len + dx2/len2 + dx3/len3 + dx4/len4);
        dplxi += (tmp1[0] + tmp1[1] + tmp1[2] + tmp1[3]);
//         dplyi += (dy/len + dy2/len2 + dy3/len3 + dy4/len4);
        dplyi += (tmp2[0] + tmp2[1] + tmp2[2] + tmp2[3]);
        
//         dplx[j]   -= tmp1[0];
//         dplx[j+1] -= tmp1[2];
//         dplx[j+2] -= tmp2[0];
//         dplx[j+3] -= tmp2[2];
//         dply[j]   -= tmp1[1];
//         dply[j+1] -= tmp1[3];
//         dply[j+2] -= tmp2[1];
//         dply[j+3] -= tmp2[3];
      }
      // and the last one is again separte
      for( ; j < nodes; ++j ) {
        dx = pxi - position[j*dim];
        dy = pyi - position[j*dim+1];
        len = dx*dx + dy*dy;
        
        // avoid division by zero
        if( len == 0. ) {
          dx = rand() * 1e-9;
          dy = rand() * 1e-9;
          len = dx*dx + dy*dy;
        }
        
        dplxi += dx/len;
        dplyi += dy/len;
        dplx[j] -= dx/len;
        dply[j] -= dy/len;
      }
      
      dplx[i] = dplxi;
      dply[i] = dplyi;
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
  aligned_free( tmp1 );
  aligned_free( tmp2 );
  
  gettimeofday(&end, NULL);
  double time_elapsed = gettime(start,end);
  printf("Layouting this graph took %1.4f sec\n", time_elapsed);
}
*/
