#ifndef __PARALLEL_REBINNING_CB_CURVE_H__
#define __PARALLEL_REBINNING_CB_CURVE_H__

extern "C"
void ParallelRebinningCBCurve_GPU(float* outputProj, 
    const float* Proj, const int YL, const int ZL, const int ViewN,
	const float YLC, const float dYA, const float DeltaTheta,
	const float PLC, const float DeltaT, const float DeltaFai,	const float SO);

#endif