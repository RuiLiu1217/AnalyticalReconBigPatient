function PProj = Parallel_Rebinning_CB_Curve(Proj,ScanGeom, useGPU)
% Rebin the spiral cone-beam projections into cone-parallel geometry
%ScanGeom = struct( 'ScDet',  [SO DO YL ZL DecAngle DecHeigh YLC ZLC h],  ...
%                   'Proj',   [BetaS BetaE ViewN N_2pi], ...
%                   'Obj',    [ObjR XN YN ZN XC YC ZC]...
%                   'Rec',    [delta HSCoef k1]);

SO = ScanGeom.ScDet(1);

YL = ScanGeom.ScDet(3);
ZL = ScanGeom.ScDet(4);
dYA= ScanGeom.ScDet(5)/YL;
YLC= ScanGeom.ScDet(7);
ViewN   = ScanGeom.Proj(3);
DeltaFai= (ScanGeom.Proj(2)-ScanGeom.Proj(1))/(ViewN-1);
DeltaTheta = DeltaFai;
DeltaT     = tan(ScanGeom.ScDet(5)*0.5)*SO*1.5*2/YL;
PLC = (1+YL)*0.5;

nProj = single(permute(Proj,[3 2 1]));
nPProj = ParallelRebinningCBCurve(single(nProj), int32(YL), int32(ZL), int32(ViewN), single(YLC), single(dYA), ...
    single(DeltaTheta), single(PLC), single(DeltaT), single(DeltaFai),single(SO),int32(useGPU));
PProj = double(permute(nPProj,[3 2 1]));

end
