% Demo of the reconstruction in CPU / GPU
GEProjectionDataName = '/home/liurui/Documents/MATLAB/BigPatient/projectionData/2117-7-1.prep.dat';
GEProjHeadDataName ='/home/liurui/Documents/MATLAB/BigPatient/projectionData/2117-7-1-hdr.dat';
ObjR = 250; % Unit: mm
Scale = 128; % size of the volume
SegSize = 800; % Each time use 800 projections
useGPU = 1; % Default use GPU

[RecImg, reconTime] = analyticalRecon(GEProjectionDataName, GEProjHeadDataName, ObjR, Scale, SegSize, useGPU);
