%------------------------------------------------------------------------- 
% Wake Forest University Health Sciences
% Univeristy of Massachusetts Lowell
% Date: 2016-08-04
% Routine: compileCPP.m
% Author
%	Rui Liu
% Email: liurui1217@gmail.com
%-------------------------------------------------------------------------
arch = computer;

if strcmp(arch, 'GLNXA64') % LINUX
    NVCC = '/usr/local/cuda/bin/nvcc';
    MATLAB_INC = [matlabroot, '/extern/include'];
    CUDA_INC = '/usr/local/cuda/include';
    CUDA_SAMP_INC = '/usr/local/cuda/samples/common/inc';
elseif  strcmp(arch,'PCWIN64') % WINDOWS
    CUDA_PATH = getenv('CUDA_PATH');
    NVCC = ['"',CUDA_PATH,'\bin\nvcc.exe"'];
    MATLAB_INC = [getenv('MATLAB_PATH'),'\extern\include'];
    CUDA_INC = [getenv('CUDA_PATH'),'\include'];
    CUDA_LIB = ['"',getenv('CUDA_PATH'),'\lib64"'];
    CUDA_SAMP_INC = [getenv('NVCUDASAMPLES_ROOT'),'\common\inc'];
elseif strcmp(arch, 'MACI64') % MAC (IS NOT IMPLEMENTED)

end
ResultName1 = 'ParallelFDKHelical3DWeightingGPU.o';
ResultName2 = 'filtering.o';
ResultName3 = 'ParallelRebinningCBCurveGPU.o';

SourceFile1 = 'ParallelFDKHelical3DWeightingGPU.cu';
SourceFile2 = 'filtering.cu';
SourceFile3 = 'ParallelRebinningCBCurveGPU.cu';

CPPFile1 = 'Parallel_FDK_Helical_3DWeighting.cpp';
CPPFile2 = 'reWeigAdFiltr_GPU.cpp';
CPPFile3 = 'ParallelRebinningCBCurve';

sys1 = [NVCC, ' -Xcompiler -O3 --use_fast_math --compile -o ', ResultName1, '  --compiler-options -fPIC -I"', MATLAB_INC, '" -I"', CUDA_INC, '" -I"', CUDA_SAMP_INC, '" ', SourceFile1];
sys2 = [NVCC, ' -Xcompiler -O3 --use_fast_math --compile -o ', ResultName2, '  --compiler-options -fPIC -I"', MATLAB_INC, '" -I"', CUDA_INC, '" -I"', CUDA_SAMP_INC, '" ', SourceFile2];
sys3 = [NVCC, ' -Xcompiler -O3 --use_fast_math --compile -o ', ResultName3, '  --compiler-options -fPIC -I"', MATLAB_INC, '" -I"', CUDA_INC, '" -I"', CUDA_SAMP_INC, '" ', SourceFile3];

% CPU/GPU based backprojection routine
system(sys1);
system(sys2);
system(sys3);
if strcmp(arch, 'GLNXA64') % LINUX
    mex -v -largeArrayDims  COMPFLAGS="$COMPFLAGS -Wall -std=c++11" -L"/usr/local/cuda/lib64" -lcudart Parallel_FDK_Helical_3DWeighting.cpp ParallelFDKHelical3DWeightingGPU.o;
    mex -v -largeArrayDims  COMPFLAGS="$COMPFLAGS -std=c++11" -L"/usr/local/cuda/lib64" -lcudart -lcufft reWeigAdFiltr_GPU.cpp filtering.o;
    mex -v -largeArrayDims  COMPFLAGS="$COMPFLAGS -std=c++11" -L"/usr/local/cuda/lib64" -lcudart -lcufft ParallelRebinningCBCurve.cpp ParallelRebinningCBCurveGPU.o;
elseif strcmp(arch, 'PCWIN64') % WINDOWS
    
    %mex -v -largeArrayDims COMPFLAGS="$COMPFLAGS -Wall -std=c++11" -L"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v8.0\lib\x64" -lcudart Parallel_FDK_Helical_3DWeighting.cpp ParallelFDKHelical3DWeightingGPU.o;
elseif strcmp(arch, 'MACI64')
    
end