#ifndef _PARALLELFDKHELICAL3DWEIGHTINGGPU_H_
#define _PARALLELFDKHELICAL3DWEIGHTINGGPU_H_

extern "C"
void ParallekFDKHelical3DWeighting_GPU_double(
		double* hctimage, const double* hfp, const double* hw, const int N_2pi,
		const int XN, const int YN, const int ZN,
		const double XNC, const double YNC, const double ZNC,
		const double dx, const double dy, const double dz,
		const int YL, const int ZL, const int PN,
		const double YLC, const double ZLC, const double dYL, const double dZL,
		const double h, const double BetaS, const double DeltaFai,
		const double N_pi, const double HSCoef, const double SO, const double SD, const double k1);

extern "C"
void ParallekFDKHelical3DWeighting_GPU_float(
		float* hctimage, const float* hfp, const float* hw, const int N_2pi,
		const int XN, const int YN, const int ZN,
		const float XNC, const float YNC, const float ZNC,
		const float dx, const float dy, const float dz,
		const int YL, const int ZL, const int PN,
		const float YLC, const float ZLC, const float dYL, const float dZL,
		const float h, const float BetaS, const float DeltaFai,
		const float N_pi, const float HSCoef, const float SO, const float SD, const float k1);
        
#endif
