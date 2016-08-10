//////////////////////////////////
//  Imaging and Informatics Lab
//  Department of Electrical and Computer Engineering
//  University of Massachusetts Lowell
//  Parallel-cone backprojection
//  May 28, 2016
//////////////////////////////////
#include <mex.h>
#include <cmath>
#include <complex>
#include <vector>
#include <cstdlib>

#include "ParallelRebinningCBCurve.h"
extern "C"
void ParallelRebinningCBCurve_GPU(float* outputProj, 
    const float* Proj, const int YL, const int ZL, const int ViewN,
	const float YLC, const float dYA, const float DeltaTheta,
	const float PLC, const float DeltaT, const float DeltaFai,	const float SO);

// Let's make the order of the projection data ZL, YL, ViewN
template<typename T>
void ParallelRebinningCBCurve_CPU(
		T* outputProj,
		const T* Proj,
		const int YL, const int ZL, const int ViewN,
		const T YLC, const T dYA,
		const T DeltaTheta,
		const T PLC, const T DeltaT, const T DeltaFai,
		const T SO)
{
	for(int i = 0; i != ViewN; ++i)
	{
		T Theta = i * DeltaTheta; // The View for the parallel projection
		for(int j = 0; j != YL; ++j)
		{
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
			else if(FI > ViewN - 1)
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
			else if(UI > YL - 1)
			{
				IndexYU = YL - 1;
				IndexYB = YL - 1;
			}
			else
			{
				IndexYU = UI;
				IndexYB = UI - 1;
			}

			for(int k = 0; k != ZL; ++k)
			{
				outputProj[(i * YL + j) * ZL + k] =
						coeXB * coeYB * Proj[(IndexXB * YL + IndexYB) * ZL + k] +
						coeXU * coeYB * Proj[(IndexXU * YL + IndexYB) * ZL + k] +
						coeXB * coeYU * Proj[(IndexXB * YL + IndexYU) * ZL + k] +
						coeXU * coeYU * Proj[(IndexXU * YL + IndexYU) * ZL + k];
			}
		}
	}
}



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // The projection is ordered in ZL, YL, ViewN order
    float* Proj = (float*)mxGetPr(prhs[0]);
    int YL = *((int*)(mxGetPr(prhs[1])));
    int ZL = *((int*)(mxGetPr(prhs[2])));
    int ViewN = *((int*)(mxGetPr(prhs[3])));
    float YLC = *((float*)(mxGetPr(prhs[4])));
    float dYA = *((float*)(mxGetPr(prhs[5])));
    float DeltaTheta = *((float*)(mxGetPr(prhs[6])));
    float PLC = *((float*)(mxGetPr(prhs[7])));
    float DeltaT = *((float*)(mxGetPr(prhs[8])));
    float DeltaFai = *((float*)(mxGetPr(prhs[9])));
    float SO = *((float*)(mxGetPr(prhs[10])));
    int useGPU = *((int*)(mxGetPr(prhs[11])));
    const mwSize dims[]={ZL,YL,ViewN};
    plhs[0] = mxCreateNumericArray(3,dims,mxSINGLE_CLASS,mxREAL);
    float* outputProj = (float*)mxGetPr(plhs[0]);
    if (useGPU==1)
    {
        ParallelRebinningCBCurve_GPU(outputProj, Proj, YL, ZL,  ViewN,
		YLC - 1.0f, dYA, DeltaTheta, PLC - 1.0f, DeltaT, DeltaFai, SO);
    }
    else
    {
             ParallelRebinningCBCurve_CPU<float>(outputProj, Proj, YL, ZL,  ViewN,
		YLC - 1.0f, dYA, DeltaTheta, PLC - 1.0f, DeltaT, DeltaFai, SO);
    }     
}

