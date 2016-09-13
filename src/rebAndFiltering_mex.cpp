/**
 * Wake Forest University Health Science
 * University of Massachusetts Lowell
 *
 * Organization:
 * 	Wake Forest University
 *
 * 	rebAndFiltering_mex.cpp
 * 	Matlab mex routine for the GPU based rebinning, reweighting
 * 	and filtering for analytical helical CT
 * 	reconstruction
 *
 * 	Author: Rui Liu
 * 	Email: liurui1217@gmail.com
 * 	Date: 2016-08-04
 *
 * 	Version 1.0
 */

#include "mex.h"
#include "matrix.h"

#include "filtering.h"

#include "ParallelRebinningCBCurve.h"
#include <iostream>
extern "C"
void filtering(float* hfpwd,
		const float* hProj,
		const int YL, const int ZL, const int ViewN,
		const float PLC, const float ZLC,
		const float dYL, const float dZL,
		const float SO);



extern "C"
void ParallelRebinningCBCurve_GPU(float* outputProj,
    const float* Proj, const int YL, const int ZL, const int ViewN,
	const float YLC, const float dYA, const float DeltaTheta,
	const float PLC, const float DeltaT, const float DeltaFai,	const float SO);

// The function routine should be called
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
	// For parallel rebinning
	float* Proj = (float*)mxGetPr(prhs[0]);
	float SO = *((float*)mxGetPr(prhs[1]));
	float DO = *((float*)mxGetPr(prhs[2]));
	int YL = *((int*)mxGetPr(prhs[3]));
	int ZL = *((int*)mxGetPr(prhs[4]));
	float dYA = *((float*)mxGetPr(prhs[5]));
	float YLC = *((float*)mxGetPr(prhs[6]));
	int ViewN = *((int*)mxGetPr(prhs[7]));
	float BetaS = *((float*)mxGetPr(prhs[8]));
	float BetaE = *((float*)mxGetPr(prhs[9]));

	float DeltaFai = (BetaE - BetaS) / (ViewN - 1);
	float DeltaTheta = DeltaFai;
	float DeltaT = tan(*((float*)mxGetPr(prhs[10])) * 0.5) * SO * 1.5 * 2.0 / YL;
	float PLC = *((float*)mxGetPr(prhs[11]));
	float* outputProj = new float[YL * ZL * ViewN];

	ParallelRebinningCBCurve_GPU(outputProj, Proj, YL, ZL, ViewN,
		YLC - 1.0f, dYA, DeltaTheta, PLC - 1.0f, DeltaT, DeltaFai, SO);


    //SO = ScanGeom.ScDet(1);
    //DO = ScanGeom.ScDet(2);
    //YL = ScanGeom.ScDet(3);
    //ZL = ScanGeom.ScDet(4);
    //dYA= ScanGeom.ScDet(5)/YL;
    //YLC= ScanGeom.ScDet(7);
    //ViewN   = ScanGeom.Proj(3);
    //DeltaFai= (ScanGeom.Proj(2)-ScanGeom.Proj(1))/(ViewN-1);
    //DeltaTheta = DeltaFai;
    //DeltaT     = tan(ScanGeom.ScDet(5)*0.5)*SO*1.5*2/YL;
    //PLC = (1+YL)*0.5;
    
    // Change the storing order
    
    // Parallel rebinning
    
    // Change the storing order
    
    // PLC, ZLC, dYL, dZL, SO are required for reweighting and filtering
    
    
    
	// YL, ZL, ViewN should be calculated from Proj
	if(nlhs != 1 || nrhs != 6)
	{
		std::cerr<<"Wrong number of parameters\n";
		std::cerr<<"This function requires 6 input parameter and provides one return value\n";
		exit(-1);
	}
	const mwSize* siz = mxGetDimensions(prhs[0]);
	const int YL = siz[0];
	const int ZL = siz[1];
	const int ViewN = siz[2];
	// Generate a new array
	plhs[0] = mxCreateNumericArray(3,siz,mxSINGLE_CLASS,mxREAL);
	float* fpwd = (float*)mxGetPr(plhs[0]);
	const float* Proj = (float*)mxGetPr(prhs[0]);
	const float PLC = *((float*)mxGetData(prhs[1]));
	const float ZLC = *((float*)mxGetData(prhs[2]));
	const float dYL = *((float*)mxGetData(prhs[3]));
	const float dZL = *((float*)mxGetData(prhs[4]));
	const float SO = *((float*)mxGetData(prhs[5]));

	filtering(fpwd, Proj, YL, ZL, ViewN,
			PLC, ZLC, dYL, dZL, SO);
    
}














