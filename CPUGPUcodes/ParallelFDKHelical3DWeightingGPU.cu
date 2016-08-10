
#include <vector>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sequence.h>
#include <cufft.h>
// includes, project
#include <cuda_runtime.h>

#include <helper_functions.h>
#include <helper_cuda.h>
#include <iostream>

template<typename T>
__global__ void ParallekFDKHelical3DWeightingKer(
		T* ctimage,
		const T* fp,
		const T* w,
		const int XN, const int YN, const int ZN,
		const T XNC, const T YNC, const T ZNC,
		const T dx, const T dy, const T dz,
		const int YL, const int ZL, const int PN,
		const T YLC, const T ZLC,
		const T dYL, const T dZL,
		const T h, const T BetaS, const T DeltaFai,
		const T N_pi, const T HSCoef, const T SO, const T SD, const T k1)
{
	int zi = threadIdx.x + blockIdx.x * blockDim.x;
	int xi = threadIdx.y + blockIdx.y * blockDim.y;
	int yi = threadIdx.z + blockIdx.z * blockDim.z;
	if(zi < ZN && xi < XN && yi < YN)
	{
		const T x = (xi - XNC) * dx;
		const T y = (yi - YNC) * dy;
		const T z = (zi - ZNC) * dz;

		// Calculate Beta0 according to the position
		const T Beta0 = 3.14159265358979323846264*2.0 * z / h;
		const int s0 = ceil((Beta0 - BetaS) / DeltaFai - 0.5);
		int s1 = s0 - ceil(N_pi * HSCoef);
		int s2 = s0 + ceil(N_pi * HSCoef);
		T res = 0;
		if((s1 < PN) && (s2 > 0))
		{
			s1 = (s1 < 0)?0:s1;
			s2 = (s2 > PN-1)?PN-1:s2;
			for(int ProjInd = 0; ProjInd <= s2; ++ProjInd)
			{
				const T View = BetaS + ProjInd * DeltaFai;
				const int d1 = N_pi - (s0 - ProjInd);
				const int d2 = (ProjInd < s0)?(d1 + N_pi) : (d1 - N_pi);
				const T UU = -x * cos(View) - y * sin(View);
				const T Yr = -x * sin(View) + y * cos(View);
				const T temp1 = sqrt(SO * SO - Yr * Yr);
				const T temp2 = (z - h * View / (2.0 * 3.14159265358979323846264));
				const T Zr = temp2*(SO*SD)/(UU * temp1+SO*SO-Yr*Yr);
				const T U1 = Yr/dYL+YLC;
				const int U  = ceil(U1);
				const T V1 = Zr/dZL+ZLC;
				const int V  = ceil(V1);
				const T Dey = U-U1;
				const T Dez = V-V1;
				if ((U>0)&&(U<YL)&&(V>0)&&(V<ZL))
				{
					const T touying =
							Dey *          Dez * fp[(ProjInd * ZL + (V-1)) * YL + (U-1)] +
							Dey * (1.0 - Dez) * fp[(ProjInd * ZL + (V)) * YL + (U-1)] +
							(1.0 - Dey) * Dez * fp[(ProjInd * ZL + (V-1)) * YL + U] +
							(1.0 - Dey) * (1.0 - Dez) * fp[(ProjInd * ZL + V) * YL + U];
					const T weight1 = w[d1];
					const T weight2 = w[d2];

					const T Gama   = fabs( temp2 / ( temp1 + UU));
					const T Gama_C = (ProjInd < s0) ? fabs((temp2 - 0.5) / (temp1 - UU)) : fabs((temp2 + 0.5) / (temp1 - UU));
					const T m1 = pow(Gama,  k1);    //m1     = std::real(std::pow(Gama,k1*h));
					const T m2 = pow(Gama_C, k1);  //m2     = std::real(std::pow(Gama_C,k1*h));
					const T weight = (weight1*m2) / (weight2*m1+weight1*m2);

					res += weight*touying*DeltaFai;
				}// end if linear interpolation

			}
		}
		ctimage[(yi*XN+xi)*ZN+zi] = res;
	}

}

template<typename T>
void ParallekFDKHelical3DWeighting_GPU(
		T* hctimage,
		const T* hfp,
		const T* hw, const int N_2pi,
		const int XN, const int YN, const int ZN,
		const T XNC, const T YNC, const T ZNC,
		const T dx, const T dy, const T dz,
		const int YL, const int ZL, const int PN,
		const T YLC, const T ZLC,
		const T dYL, const T dZL,
		const T h, const T BetaS, const T DeltaFai,
		const T N_pi, const T HSCoef, const T SO, const T SD, const T k1)
{
	dim3 blk(64,16,1);
	dim3 gid(
			(ZN + blk.x - 1) / blk.x,
			(XN + blk.y - 1) / blk.y,
			(YN + blk.z - 1) / blk.z);

	thrust::device_vector<T> ctimage(XN * YN * ZN, 0);
	thrust::device_vector<T> fp(hfp, hfp + YL * ZL * PN);
	thrust::device_vector<T> w(hw, hw + N_2pi);
	ParallekFDKHelical3DWeightingKer<T><<<gid,blk>>>(
			thrust::raw_pointer_cast(&ctimage[0]),
			thrust::raw_pointer_cast(&fp[0]),
			thrust::raw_pointer_cast(&w[0]), XN, YN, ZN, XNC, YNC, ZNC,
			dx, dy, dz, YL, ZL, PN, YLC, ZLC, dYL, dZL, h, BetaS, DeltaFai, N_pi,
			HSCoef, SO, SD, k1);
	thrust::copy(ctimage.begin(),ctimage.end(),hctimage);
    ctimage.clear();
    fp.clear();
    w.clear();
}


extern "C"
void ParallekFDKHelical3DWeighting_GPU_double(
		double* hctimage, const double* hfp, const double* hw, const int N_2pi,
		const int XN, const int YN, const int ZN,
		const double XNC, const double YNC, const double ZNC,
		const double dx, const double dy, const double dz,
		const int YL, const int ZL, const int PN,
		const double YLC, const double ZLC, const double dYL, const double dZL,
		const double h, const double BetaS, const double DeltaFai,
		const double N_pi, const double HSCoef, const double SO, const double SD, const double k1)
{
	ParallekFDKHelical3DWeighting_GPU<double>(hctimage, hfp, hw, N_2pi,
			XN, YN, ZN, XNC, YNC, ZNC, dx, dy, dz,
			YL, ZL, PN, YLC, ZLC, dYL, dZL, h, BetaS, DeltaFai,
			N_pi, HSCoef, SO, SD, k1);
}



extern "C"
void ParallekFDKHelical3DWeighting_GPU_float(
		float* hctimage, const float* hfp, const float* hw, const int N_2pi,
		const int XN, const int YN, const int ZN,
		const float XNC, const float YNC, const float ZNC,
		const float dx, const float dy, const float dz,
		const int YL, const int ZL, const int PN,
		const float YLC, const float ZLC, const float dYL, const float dZL,
		const float h, const float BetaS, const float DeltaFai,
		const float N_pi, const float HSCoef, const float SO, const float SD, const float k1)
{
	ParallekFDKHelical3DWeighting_GPU<float>(hctimage, hfp, hw, N_2pi,
			XN, YN, ZN, XNC, YNC, ZNC, dx, dy, dz,
			YL, ZL, PN, YLC, ZLC, dYL, dZL, h, BetaS, DeltaFai,
			N_pi, HSCoef, SO, SD, k1);
}
