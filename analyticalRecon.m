% Author: Rui Liu
% Email: liurui1217@gmail.com
% Organization: Wake Forest Health Sciences & University of Massachusetts
% Lowell
% Routine: The CPU/GPU routine for big patient reconstruction
% NOTE: The current script only support the GE projection data. 
function [RecImg, reconTime] = analyticalRecon(GEProjectionDataName, GEProjHeadDataName, ObjR, volumeSize, SegLength, useGPU)
%% Read the projection data
close all;
clc;

%% Recon the volume
objR = ObjR;
Scale= volumeSize;
xn = Scale;
yn = Scale;
zn = Scale;%ceil(Scale*0.5);
mWtInd = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CTGE = GetHeadofGeDataSet(GEProjHeadDataName);
SegOverLap = atan(CTGE.DetCellPitch/CTGE.Source2DetDis/2)*CTGE.NumberOfDetChannelPerRow;
SegOverLap = ceil(SegOverLap*CTGE.NumberOfViewPerRotation/(2*pi))+1;
SegNum = ceil((double(CTGE.TotalViewNumber)-SegLength)/(SegLength-2*SegOverLap))+1;
ctimage = zeros(xn,yn,zn,SegNum);
disp('Totally we divide the number of Segments with ');
disp(SegNum);
tic;
for SegInd = 1:SegNum
  [Proj, DLength] = readBigPatientData(GEProjectionDataName, CTGE, SegLength, SegOverLap, SegInd);
  [res, recontime] = analyticalReconKer(Proj, CTGE, SegLength, SegOverLap, SegInd, DLength, objR, xn, yn, zn, mWtInd, useGPU);
  ctimage(:,:,:,SegInd) = res;
end;
RecImg = sum(ctimage,4);

reconTime = toc;
imdisp(RecImg);

end