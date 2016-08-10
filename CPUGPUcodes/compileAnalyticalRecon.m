%------------------------------------------------------------------------- 
% Wake Forest University Health Sciences
% Univeristy of Massachusetts Lowell
% Date: 2016-08-04
% Routine: compileCPP.m
% Author
%	Rui Liu
% Email: liurui1217@gmail.com
%-------------------------------------------------------------------------

% CPU/GPU based backprojection routine
system( '/usr/local/cuda/bin/nvcc -Xcompiler -O3 --use_fast_math --compile -o ParallelFDKHelical3DWeightingGPU.o  --compiler-options -fPIC  -I"/usr/local/MATLAB/R2015b/extern/include " -I/usr/local/cuda/include -I/usr/local/cuda/samples/common/inc "ParallelFDKHelical3DWeightingGPU.cu" ' );
mex -v -largeArrayDims  COMPFLAGS="$COMPFLAGS -std=c++11" -L/usr/local/cuda/lib64 -L/usr/local/cuda/lib64 -lcudart Parallel_FDK_Helical_3DWeighting.cpp ParallelFDKHelical3DWeightingGPU.o 

% GPU based filtering and reweighting
system( '/usr/local/cuda/bin/nvcc -lcufft -Xcompiler -O3 --use_fast_math --compile -o filtering.o  --compiler-options -fPIC  -I"/usr/local/MATLAB/R2015b/extern/include " -I/usr/local/cuda/include -I/usr/local/cuda/samples/common/inc "filtering.cu" ' );
mex -v -largeArrayDims  COMPFLAGS="$COMPFLAGS -std=c++11" -L/usr/local/cuda/lib64 -lcudart -lcufft reWeigAdFiltr_GPU.cpp filtering.o 

% GPU / CPU based rebinning
system( '/usr/local/cuda/bin/nvcc -lcufft -Xcompiler -O3 --use_fast_math --compile -o ParallelRebinningCBCurveGPU.o  --compiler-options -fPIC  -I"/usr/local/MATLAB/R2015b/extern/include " -I/usr/local/cuda/include -I/usr/local/cuda/samples/common/inc "ParallelRebinningCBCurveGPU.cu" ' );
mex -v -largeArrayDims  COMPFLAGS="$COMPFLAGS -std=c++11" -L/usr/local/cuda/lib64 -lcudart -lcufft ParallelRebinningCBCurve.cpp ParallelRebinningCBCurveGPU.o 
