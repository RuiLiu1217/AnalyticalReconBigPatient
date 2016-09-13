% -------------------------------------------------------------------------
% Author: Rui Liu
% Date: May. 19, 2016
% Wake Forest Health Sciences & University of Massachusetts Lowell
% The routine to read the GE data in "dat" form and save it to a raw form
% that can be reconstructed by GE_DD3 function in CUDA. 
% Input:
%       ProjName        :          Name of projection data in dat format
%       ProjHeadName    :          Head file name of projection data in
%       dat format
%       GEProjFileName  :          The file name that to save for GE_DD3
%       projection routine.
% Output: 
%       Proj 
function [Proj,DLength] = readBigPatientData(ProjName,CTGE, SegLength,SegOverLap, SegInd)
SegNum = ceil((double(CTGE.TotalViewNumber)-SegLength)/(SegLength-2*SegOverLap))+1;
if (SegInd<SegNum)
    DLength = SegLength;
else
    DLength = CTGE.TotalViewNumber-(SegLength-2*SegOverLap)*(SegInd-1);
end;
FID = fopen(ProjName);
SkipNum = (SegLength-2*SegOverLap) * CTGE.NumberOfDetChannelPerRow * CTGE.NumberOfDetectorRow * (SegInd-1);
Status =fseek(FID,SkipNum*4,'bof');
if (Status>=0)
  Proj = fread(FID, CTGE.NumberOfDetChannelPerRow*CTGE.NumberOfDetectorRow*DLength, 'float');
  Proj = reshape(Proj, CTGE.NumberOfDetChannelPerRow, CTGE.NumberOfDetectorRow, DLength);
  Proj = permute(Proj,[3 1 2]);
  %Proj = Proj(:,end:-1:1,:);
else
  Proj = -1;
end;

fclose(FID);

