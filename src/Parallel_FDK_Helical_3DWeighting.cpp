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
//#define DOUBLE_PRECISION 
#define TYPE double
const static TYPE pi = 3.14159265358979;

#include "ParallelFDKHelical3DWeightingGPU.h"

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

void ParallelFDKHelical3DWeightingCPU(double* ctimage, const double* fp, const double* w, const int N_2pi,
		const int XN, const int YN, const int ZN,
		const double XNC, const double YNC, const double ZNC,
		const double dx, const double dy, const double dz,
		const int YL, const int ZL, const int PN,
		const double YLC, const double ZLC, const double dYL, const double dZL,
		const double h, const double BetaS, const double DeltaFai,
		const double N_pi, const double HSCoef, const double SO, const double SD, const double k1)
{   
    const double invPI = 1.0 / (2.0 * 3.141592653589793);
    for(int yi=0;yi<YN;yi++)
    {
        const double y = (yi-YNC)*dy;
        for(int xi=0;xi<XN;xi++)
        {
            const double x  = (xi-XNC)*dx;
            for(int zi = 0; zi < ZN; ++zi)
            {
                    ///compute the projection position for every grid on the image plane
                    const double z = (zi-ZNC) * dz;
                    const double Beta0 = 2.0 * pi * z / h;
                    const int s0 = ceil((Beta0-BetaS) / DeltaFai-0.5);
                    int s1 = s0-ceil(N_pi*HSCoef);       //s1 = ceil((Beta0 + N_Circle * pi - pi ) /DeltaFai);
                    int s2 = s0+ceil(N_pi*HSCoef);       //s2 = ceil((Beta0 + N_Circle * pi + pi ) /DeltaFai);
                    double res = 0; // used to accumulate the results
                    if((s1 < PN) && (s2 > 0))
                    {
                        s1 = (s1 < 0)?0:s1;
                        s2 = (s2 > PN - 1)?PN-1:s2;
                        for(int ProjInd = s1; ProjInd <= s2; ++ProjInd)
                        {
                                const double View = BetaS + ProjInd * DeltaFai;
                                const int d1   = N_pi-(s0-ProjInd); //d1 = ProjInd;
                                const int d2 = (ProjInd < s0)?(d1 + N_pi) : (d1 - N_pi);
                                const double UU = -x*cos(View)-y*sin(View);
                                const double Yr = -x*sin(View)+y*cos(View);
                                const double temp1 = sqrt(SO * SO - Yr * Yr);
                                const double temp2 = (z-h*View * invPI);
                                const double Zr = temp2*(SO*SD)/(UU * temp1+SO * SO - Yr * Yr);
                                const double U1 = Yr/dYL+YLC;
                                const int U  = ceil(U1);
                                const double V1 = Zr/dZL+ZLC;
                                const int V  = ceil(V1);
                                const double Dey = U-U1;
                                const double Dez = V-V1;
                                //Linear interploate
                                if ((U>0)&&(U<YL)&&(V>0)&&(V<ZL))
                                {
                                    const double touying =
                                            Dey *          Dez * fp[(ProjInd * ZL + (V-1)) * YL + (U-1)] +
                                            Dey * (1.0 - Dez) * fp[(ProjInd * ZL + (V)) * YL + (U-1)] +
                                            (1.0 - Dey) * Dez * fp[(ProjInd * ZL + (V-1)) * YL + U] +
                                            (1.0 - Dey) * (1.0 - Dez) * fp[(ProjInd * ZL + V) * YL + U];
                                    const double weight1 = w[d1];
                                    const double weight2 = w[d2];

                                    const double Gama   = fabs( temp2 / ( temp1 + UU));
                                    const double Gama_C = (ProjInd < s0) ? fabs((temp2 - 0.5) / (temp1 - UU)) : fabs((temp2 + 0.5) / (temp1 - UU));
                                    const double m1 = pow(Gama,  k1);    //m1     = std::real(std::pow(Gama,k1*h));
                                    const double m2 = pow(Gama_C, k1);  //m2     = std::real(std::pow(Gama_C,k1*h));
                                    const double weight = (weight1*m2) / (weight2*m1+weight1*m2);

                                    res += weight*touying*DeltaFai;
                                }// end if linear interpolation
                        }// end for projection

                    }// end if range
                    ctimage[(yi*XN+xi)*ZN+zi] = res;
            } // end zi
        }	// xi
    }// yi
}



extern "C"
void ParallelFDKHelical3DWeighting(
		double* ctimage, // Image to be reconstructed
		const double* fp,
		const double SO,
		const double DO,
		const int YL,
		const int ZL,
		const double DecWidth,
		const double DecHeigh,
		const double YLC,
		const double ZLC,
		const double h,
		const double BetaS,
		const int PN,
		const int N_2pi,
		const double ObjR,
		const int XN,
		const int YN,
		const int ZN,
		const double XNC,
		const double YNC,
		const double ZNC,
		const int delta,
		const double HSCoef,
		const double k1,
        const int useGPU)
{
        const double PI = 3.141592653589793;
        const double N_pi = N_2pi / 2.0;
        const double dx = 2.0 * ObjR / XN;
        const double dy = 2.0 * ObjR / YN;
        const double dz = 2.0 * ObjR / ZN;
        const double dYL = DecWidth / YL;
        const double dZL = DecHeigh / ZL;
        const double DeltaFai = 2.0 * PI / N_2pi;
        const double inv2PI = 1.0 / (2.0 * PI);
        const double SD = SO + DO;
        const double SO_square = SO * SO;
        std::vector<double> w(N_2pi,0);
        const int L = 2.0 * ceil(N_pi * HSCoef) + 1.0;
        const int Shift = N_pi - ceil(N_pi * HSCoef);
    
        ////Producing the weighting function

        for (int k=0;k<L;k++)
        {
            if (0 <= k && k<delta)
                w[k+Shift]= pow(cos((pi/2)*(delta-k-0.5)/delta),2);
            else if(L-delta<=k && k < L)
                w[k+Shift]= pow(cos((pi/2)*(k-(L-delta)+0.5)/delta),2);
            else
                w[k+Shift] = 1;
        }

        
        if(useGPU==1)
        {

            // use single floating data to implement the backprojection
             std::vector<float> fw(w.begin(),w.end());
             float* fctimage = new float[XN * YN * ZN];

             float* ffp = new float[YL * ZL * PN];
             std::copy(fp, fp + YL * ZL * PN, ffp);

             ParallekFDKHelical3DWeighting_GPU_float(fctimage,ffp, &fw[0],N_2pi,
                     XN, YN, ZN, XNC, YNC, ZNC, dx, dy, dz,
                      YL, ZL, PN, YLC, ZLC, dYL, dZL, h, BetaS, DeltaFai,
                      N_pi, HSCoef, SO, SD, k1);
             std::copy(fctimage, fctimage + XN * YN * ZN, ctimage);
             delete[] fctimage;
             delete[] ffp;
        }
        else
        {
            ParallelFDKHelical3DWeightingCPU(ctimage, fp, &w[0], N_2pi,
                    XN, YN, ZN, XNC, YNC, ZNC, dx, dy, dz,
                     YL, ZL, PN, YLC, ZLC, dYL, dZL, h, BetaS, DeltaFai,
                     N_pi, HSCoef, SO, SD, k1);
        }
        
        //////end of the main code
}


static const TYPE *geom[4];
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    TYPE ObjR,SO,DO,YLC,ZLC,XNC,YNC,ZNC,DecWidth,DecHeigh,h,dYL,dZL,dx,dy,dz,DeltaFai,BetaS,k1,HSCoef;
    int  YL,ZL,XN,YN,ZN,N_2pi,delta,N_pi,PN,i,j,k;
    TYPE *ctimage,*w;
    const TYPE *fp;
    /* Check for proper number of arguments */
    if (nrhs != 4) {
        mexErrMsgTxt("Backward projection needs 4 inputs.");
    }
    if (nlhs != 0) {
        mexErrMsgTxt("Backward projection does not need any output.");
    }
    //Parallel_FDK_Helical_3DWeighting(ScanGeom,fpwd,ctimage);
    //ScanGeom = struct( 'ScDet',  [SO DO YL ZL DecAngle DecHeigh YLC ZLC h],  ...
    //               'Proj',   [BetaS BetaE ViewN N_2pi], ...
    //               'Obj',    [ObjR XN YN ZN XC YC ZC], ...
    //               'Rec',    [delta HSCoef k1]);
    geom[0] = mxGetPr(mxGetFieldByNumber(prhs[0], 0, 0));
    SO            = TYPE(geom[0][0]);
    DO            =    TYPE(geom[0][1]);
    YL            = int(TYPE(geom[0][2])+0.5);
    ZL            = int(TYPE(geom[0][3])+0.5);
    DecWidth      = TYPE(geom[0][4]);
    DecHeigh      = TYPE(geom[0][5]);
    YLC           =     (YL-1)*0.5;//geom[0][6]-1;
    ZLC           =   TYPE(geom[0][7])-1.0;
    h             =     TYPE(geom[0][8])*DecHeigh;
    
    geom[1] = mxGetPr(mxGetFieldByNumber(prhs[0], 0, 1));
    BetaS         =     TYPE(geom[1][0]);
    PN            =     int(geom[1][2]);
    N_2pi         = int(geom[1][3]);
    
    geom[2] = mxGetPr(mxGetFieldByNumber(prhs[0], 0, 2));
    ObjR          =    TYPE(geom[2][0]);
    XN            = int(TYPE(geom[2][1])+0.5);
    YN            = int(TYPE(geom[2][2])+0.5);
    ZN            = int(TYPE(geom[2][3])+0.5);
    XNC           =     TYPE(geom[2][4])-1.0;
    YNC           =     TYPE(geom[2][5])-1.0;
    ZNC           =     TYPE(geom[2][6])-1.0;
    
    geom[3] = mxGetPr(mxGetFieldByNumber(prhs[0], 0, 3));
    delta         = int(TYPE(geom[3][0])+0.5);
    HSCoef        =     TYPE(geom[3][1]);// It is the value of Half_Scan/2*pi;
    k1            =     TYPE(geom[3][2])*TYPE(geom[0][8]);
    
    fp            =     mxGetPr(prhs[1]);
    ctimage       =     mxGetPr(prhs[2]);
    int useGPU = *((int*)mxGetPr(prhs[3]));
    
    
    N_pi = N_2pi/2.0;
    dx = 2.0*ObjR/XN;
    dy = 2.0*ObjR/YN;
    dz = 2.0*ObjR/ZN;
    dYL = DecWidth/YL;
    dZL = DecHeigh/ZL;
    DeltaFai = 2.0*pi/N_2pi;
    
     ParallelFDKHelical3DWeighting(ctimage, fp, SO, DO, YL, ZL,
		DecWidth, DecHeigh, YLC, ZLC, h, BetaS, PN, N_2pi,
		ObjR, XN, YN, ZN, XNC, YNC, ZNC, delta, HSCoef, k1,useGPU);
    
}

