function [RecImg, RecTime] = analyticalRecon(GEProjectionDataName, GEProjHeadDataName, SegLength, objR, Scale, useGPU)

% %% Read the projection data
% clear all;
close all;
clc;

xn = Scale;
yn = Scale;
zn = Scale;%ceil(Scale*0.5);
mWtInd = 3; % use 3D weighting

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CTGE = GetHeadofGeDataSet(GEProjHeadDataName);
SegOverLap = atan(CTGE.DetCellPitch/CTGE.Source2DetDis/2)*CTGE.NumberOfDetChannelPerRow;
SegOverLap = ceil(SegOverLap*CTGE.NumberOfViewPerRotation/(2*pi))+1;
SegNum = ceil((double(CTGE.TotalViewNumber)-SegLength)/(SegLength-2*SegOverLap))+1;
RecImg = zeros(xn,yn,zn);
totaltime =0;
for SegInd = 1:SegNum
  [Proj DLength] = readBigPatientData(GEProjectionDataName, CTGE, SegLength, SegOverLap, SegInd);
  [ctimage, recontime] = analyticalReconBigPatient(Proj, CTGE, SegLength, SegOverLap, SegInd, DLength, objR, xn, yn, zn, mWtInd, useGPU);
  RecImg = RecImg+ctimage;
  totaltime = totaltime+recontime;
end;
RecTime = totaltime;

end

