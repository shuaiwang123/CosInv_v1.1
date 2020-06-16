%function MakeDtopo
% FUNCTION MAKEDTOPO
%
% This code is used to prepare the dtopo files (dtopotype=3) required by
% Geoclaw. For tsunami modelling a file dtopo file is generally used to
% specify the displacement of the topography relative to that specified in
% the topo files. This code will make dtopo file(seafloor deformation) for
% 1 m of strike slip and 1 m of dip slip motion on all subfaults. If you
% have done this you can run the python script RUNGEOCLAW.py to make system
% call of Geoclaw software.
%
% OUT
%
% NOTE
% The principle of this script is from Diego. However, some comments have
% been added by me, and this new version of code looks like more beautiful
% and more easy to use. If you have any questions or suggestions about this
% script, please do not hesitate to email me, I do like to reply you.
%
% STEP
% 0. Pre-settings
% 1. Make the seafloor pointws
% 2. Load fault model
% 3. Make GREENs for seafloor points
% 4. Prepare the dtopo file required by geoclaw
%
% By Shuai WANG @PolyU 2017-08-16
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clc;clear;
disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
disp('++++++++++++++++      Welcome To MakeDtopo         ++++++++++++++++');
disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
%==========================================================================
%================                  0                       ================
%======                       Pre-settings                           ======
%==========================================================================
disp('...... Pre-settings ......');
% addpath('Irian Coseismic Disp');
origin        = [-0.51 132.78];
code_dir      = '/media/wenbin/Student/shuaiw/CosInv_v1.1/Code';
SetPaths(code_dir);
addpath('/media/wenbin/Student/shuaiw/CosInv_v1.1/Irian Coseismic Disp');

home          = '/media/wenbin/Student/shuaiw/CosInv_v1.1';
project_name  = 'IrianTsu';
geoclaw_name  = 'geoclaw';
run_name      = 'Irian';

% Spatial range to create the seafloor grid points
lon_min       = 130;
lon_max       = 135;
lat_min       = -0.5;
lat_max       = 1;
lon_mesh_cre      = 0.0167;
lat_mesh_cre      = 0.0167;

%==========================================================================
%================                  1                       ================
%======                Make the Seafloor Points                      ======
%==========================================================================
MakeSFPoints;


%==========================================================================
%================                  2                       ================
%======                     Load Fault Model                         ======
%==========================================================================
LoadFaultModel;


%==========================================================================
%================                  3                       ================
%======                Make GREENs for Seafloor                      ======
%==========================================================================
MakeSFResponse;


%==========================================================================
%================                  4                       ================
%======          Prepare the Dtopo Required by Geoclaw               ======
%==========================================================================
Prep4Dtopo;










