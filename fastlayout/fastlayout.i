%module fastlayout

%{
    #define SWIG_FILE_WITH_INIT
    #include "fastlayout.h"
%}

%include "numpy.i"

%init %{
    import_array();
%}


%apply (float* ARGOUT_ARRAY1, int DIM1) {
(float *position, int nodes)
}
%apply (double* ARGOUT_ARRAY1, int DIM1) {
(double *position, int nodes)
}
%apply (int* IN_ARRAY2, int DIM1, int DIM2) {
(int *es, int edges, int dummy_dim)
}


%include "fastlayout.h"