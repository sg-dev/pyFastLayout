%module fastlayout

%{
    #define SWIG_FILE_WITH_INIT
    #include "fastlayout.h"
%}

%include "numpy.i"

%init %{
    import_array();
%}

%apply (foat* INPLACE_ARRAY1, int DIM1) {
(float *position, int nodes)
}

%apply (float* INPLACE_ARRAY2, int DIM1, int DIM2) {
(float *position, int nodes, int dim)
}
%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2) {
(double *position, int nodes, int dim)
}
%apply (int* IN_ARRAY2, int DIM1, int DIM2) {
(int *es, int edges, int dummy_dim)
}


%include "fastlayout.h"