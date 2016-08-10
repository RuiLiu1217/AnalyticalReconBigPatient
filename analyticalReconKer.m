%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is to demonstrate the application of analytic cone-beam reconstrution with curved detector. 
% The codes were developed by Hengyong Yu (Department of Electrical and Computer Engineering, 
% University of Massachusetts Lowell) and Rui Liu (Wake Forest Health Sciences).
% Email: liurui1217@gmail.com
% Version 2.0
% Date: Aug. 10, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ctimage, time] = analyticalReconKer(Proj, CTGE, SegLength, SegOverLap, SegInd, DLength, ...
    objR, xn, yn, zn, mWtInd, useGPU, ...
    mDetHei, mk1, mdelta, mHSCoef)

 if(~exist('mDetHei','var') || isempty(mDetHei))
     mDetHei = 1.0963;
 end

if(~exist('mk1','var') || isempty(mk1))
    mk1 = 5;
end

if(~exist('mdelta','var') || isempty(mdelta))
    mdelta = 60;
end 

if(~exist('mHSCoef','var') || isempty(mHSCoef))
    mHSCoef = 0.68;
end

Proj = Proj(:,end:-1:1,:);
CTGE.IsoChannelLocation = CTGE.NumberOfDetChannelPerRow-CTGE.IsoChannelLocation;
time = cputime;
% Rebinning (Convert conebeam to parallel beam)
% The following lines 13-23 are used to define the helical scanning geometry for a HD Geometry CT 
SO = CTGE.Source2IsoDis;          % The distance from source to object (cm)
DO = CTGE.Source2DetDis - CTGE.Source2IsoDis;   % The distance from detector to object center(cm)
YL = CTGE.NumberOfDetChannelPerRow;         % Detector cell number along the horizontal direction of detector array
YLC= CTGE.IsoChannelLocation;  % Detector center along the horizontal direction of detector array
ZL = CTGE.NumberOfDetectorRow;          % Detector cell number along the vertical direction of detector array
ZLC= (1+ZL)*0.5;  %Detector center along the vertical direction of detector array 
                  % TODO: Generalize this parameter
DetWid = CTGE.DetCellPitch; % Detector size in in-plane direciton (Set by Rui Liu)
%DetHei = mDetHei;%
DetHei = CTGE.SliceHeigh*(DO+SO)/SO; % //M by Yu //Detector size in cross-plane direction (Set by Rui Liu)

DecAngle = atan(DetWid / (SO + DO) / 2) * YL * 2; % Detector array beam angle along the horizontal direction (rad)
DecHeigh = DetHei * ZL;   % Detector array height along the vertical direction (cm)
%N_Turn   = 8;     % The number of turns for the whole helical scan 
N_2pi    = CTGE.NumberOfViewPerRotation;   % The projections/views number for each turn of scan
h        = CTGE.HelicalPitch / CTGE.NumberOfDetectorRow*SO/(DO+SO); %//M by Yu        % Helical pitch related to detector height


%The following line 25-28 are used to define the reconstruction paramters
WtInd    = mWtInd;     %1:Redundancy weighting; 2:2D weigthing; 3: 3D weighting 
k1       = mk1;     % The order to define the 3D weighting function
delta    = mdelta;    % The range to define smoothness of 2D weigthing function 
HSCoef   = mHSCoef;   % This is used to define the half-scan range              

ObjR  = objR;      % Diameter of imaginh object
%SSO   = 64;       % Define the size of the reconstructed image
XN    = xn;
YN    = yn;
ZN    = zn;

XC    = (XN+1) * 0.5;
YC    = (YN+1) * 0.5;
ZC    = 10;

ViewN = DLength;
SegNum = ceil((double(CTGE.TotalViewNumber)-SegLength)/(SegLength-2*SegOverLap))+1;
BetaS = pi-CTGE.BeginViewAngle * pi/ 180.0 +(SegInd-1)*(SegLength-2*SegOverLap)* 2 * pi / CTGE.NumberOfViewPerRotation;

BetaE = BetaS + ViewN * 2 * pi / CTGE.NumberOfViewPerRotation;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ScanGeom = struct( 'ScDet',  [SO DO YL ZL DecAngle DecHeigh YLC ZLC h],  ...
                   'Proj',   [BetaS BetaE ViewN N_2pi], ...
                   'Obj',    [ObjR XN YN ZN XC YC ZC], ...
                   'Rec',    [delta HSCoef k1]);
               
%% Purely inside GPU reconstruction codes
% NOTE: THE PRECISION IS BAD, IF WE WANT TO FURTHER IMPROVE THE SPEED PERFORMANCE,
% THE FOLLOWING TWO LINE CODES NEED TO BE DEBUGGED!!!!
% Proj = permute(Proj,[2  3  1]);
% ctimage = double(analyticalReconGPU_mex(single(Proj), single(SO),single(DO), int32(YL),int32(ZL), single(DecAngle), single(DecHeigh), single(YLC), single(ZLC), single(h), single(BetaS), single(BetaE), int32(ViewN), int32(N_2pi),single(ObjR), int32(XN), int32(YN), int32(ZN), single(XC), single(YC), single(ZC), single(delta), single(HSCoef), single(k1),int32(SegOverLap)));
    
              

% Different part end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rebinning conebeam to parallel beam 
disp('Parallel rebinning from the conebeam projections');
Proj = Parallel_Rebinning_CB_Curve(Proj,ScanGeom,useGPU); % is replaced by the codes below    

if (SegInd==1)
    Proj  = Proj(1:(ViewN-SegOverLap),:,:);
    ViewN = ViewN-SegOverLap;
    BetaE = BetaE - SegOverLap * 2 * pi / CTGE.NumberOfViewPerRotation;
elseif(SegInd==SegNum)
    Proj = Proj((SegOverLap+1):ViewN,:,:);
    ViewN = ViewN-SegOverLap;
    BetaS = BetaS + SegOverLap * 2 * pi / CTGE.NumberOfViewPerRotation;
else
    Proj = Proj((SegOverLap+1):(ViewN-SegOverLap),:,:);
    ViewN = ViewN-2*SegOverLap;
    BetaS = BetaS + SegOverLap * 2 * pi / CTGE.NumberOfViewPerRotation;
    BetaE = BetaE - SegOverLap * 2 * pi / CTGE.NumberOfViewPerRotation;
end; 

DecWidth = tan(DecAngle*0.5)*SO*2*1.5; % for rebinning
dYL   =  DecWidth/YL;
dZL   =  DecHeigh/ZL;


ScanGeom = struct( 'ScDet',  [SO DO YL ZL DecWidth DecHeigh YLC ZLC h],  ...
                   'Proj',   [BetaS BetaE ViewN N_2pi], ...
                   'Obj',    [ObjR XN YN ZN XC YC ZC], ...
                   'Rec',    [delta HSCoef k1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Preweighting the conebeam projections

disp('PreWeighting the cone beam projection');
PLC = (YL - 1)*0.5;
% currently the projection data is stored in the form 
% Proj(ViewN, Transverse, Vertical);
newProj = permute(Proj,[2 3 1]);
fpwd = double(reWeigAdFiltr_GPU(single(newProj),single(PLC), single(ZLC - 1.0), single(dYL), single(dZL), single(SO)));

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Backproject the filtered data into the 3D space
disp('BackProjection the filtered projection');
ctimage = zeros(XN,YN,ZN);
% Backprojection with GPU / CPU in 3D weighting
Parallel_FDK_Helical_3DWeighting(ScanGeom,fpwd,ctimage,int32(useGPU)); % USE GPU VERSION
time = cputime-time;

end


