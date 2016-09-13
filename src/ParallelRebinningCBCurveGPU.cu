#include "ParallelRebinningCBCurve.h"
#include <thrust/device_vector.h>
#include <thrust/copy.h>


template<typename T>
__global__ void ParallelRebinningCBCurve_GPU_ker(
		T* outputProj,
		const T* Proj,
		const int YL, const int ZL, const int ViewN,
		const T YLC, const T dYA,
		const T DeltaTheta,
		const T PLC, const T DeltaT, const T DeltaFai,
		const T SO)
{
	int k = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;
	int i = threadIdx.z + blockIdx.z * blockDim.z;
	if(i < ViewN && j < YL && k < ZL)
	{
		T Theta = i * DeltaTheta; // The View for the parallel projection
		T t = (j - PLC) * DeltaT;
		T Beta = asin(t / SO);
		T Fai = Theta + Beta;
		T a = atan(t / sqrt(SO*SO-t*t));
		T FaiIndex = (Fai / DeltaFai);
		T UIndex = a / dYA + YLC;
		int FI = ceil(FaiIndex);
		int UI = ceil(UIndex);
		T coeXB = FI - FaiIndex;
		T coeXU = 1.0 - coeXB;
		T coeYB = UI - UIndex;
		T coeYU = 1.0 - coeYB;

		int IndexXU(0);
		int IndexXB(0);
		int IndexYU(0);
		int IndexYB(0);

		if(FI <= 0)
		{
			IndexXU = 0;
			IndexXB = 0;
		}
		else if(FI >= ViewN - 1)
		{
			IndexXU = ViewN - 1;
			IndexXB = ViewN - 1;
		}
		else
		{
			IndexXU = FI;
			IndexXB = FI - 1.0;
		}

		if(UI <= 0)
		{
			IndexYU = 0;
			IndexYB = 0;
		}
		else if(UI >= YL - 1)
		{
			IndexYU = YL - 1;
			IndexYB = YL - 1;
		}
		else
		{
			IndexYU = UI;
			IndexYB = UI - 1;
		}
		outputProj[(i * YL + j) * ZL + k] =
			coeXB * coeYB * Proj[(IndexXB * YL + IndexYB) * ZL + k] +
			coeXU * coeYB * Proj[(IndexXU * YL + IndexYB) * ZL + k] +
			coeXB * coeYU * Proj[(IndexXB * YL + IndexYU) * ZL + k] +
			coeXU * coeYU * Proj[(IndexXU * YL + IndexYU) * ZL + k];
	}
}

template<typename T>
void ParallelRebinningCBCurve_GPU_template(
		T* outputProj,
		const T* Proj,
		const int YL, const int ZL, const int ViewN,
		const T YLC, const T dYA,
		const T DeltaTheta,
		const T PLC, const T DeltaT, const T DeltaFai,
		const T SO)
{
	dim3 blk(64,16,1);
	dim3 gid(
		(ZL + blk.x - 1) / blk.x,
		(YL + blk.y - 1) / blk.y,
		(ViewN + blk.z - 1) / blk.z);
	thrust::device_vector<T> dProj(Proj, Proj + YL * ZL * ViewN);
	thrust::device_vector<T> doutputProj(dProj.size(),0);
	ParallelRebinningCBCurve_GPU_ker<T><<<gid,blk>>>(
			thrust::raw_pointer_cast(&doutputProj[0]),
			thrust::raw_pointer_cast(&dProj[0]), YL, ZL, ViewN,
			YLC, dYA, DeltaTheta, PLC, DeltaT,
			DeltaFai, SO);
	thrust::copy(doutputProj.begin(),doutputProj.end(),outputProj);
	dProj.clear();
	doutputProj.clear();
}

extern "C"
void ParallelRebinningCBCurve_GPU(float* outputProj, 
    const float* Proj, const int YL, const int ZL, const int ViewN,
	const float YLC, const float dYA, const float DeltaTheta,
	const float PLC, const float DeltaT, const float DeltaFai,	const float SO)
{
    cudaSetDevice(0);
    ParallelRebinningCBCurve_GPU_template<float>(outputProj, Proj,
		YL, ZL, ViewN, YLC, dYA, DeltaTheta, PLC, DeltaT, DeltaFai, SO);
}






