% This scripts adds all the subfolders to MATLAB path, so that any function inside them can
% can be accessed from parent directory. Please run this script at the beginning of any session.

currDir = pwd;
HSfolder = strcat(currDir,'/ViscousHessSmith');
XFfolder = strcat(currDir,'/XFOIL');
OPTIfolder = strcat(currDir,'/OptimizationCode');
LAPfolder = strcat(currDir,'/Lap_Performance');
GEOMfolder = strcat(currDir,'/Geometry');
CORRfolder = strcat(currDir,'/2dto3dcorrection');
CFDfolder = strcat(currDir,'/CFD');
addpath(HSfolder,XFfolder,LAPfolder,GEOMfolder,CORRfolder,OPTIfolder,CFDfolder);
