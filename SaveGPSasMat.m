% This scripts is used to save the GPS coseismic displacement as .mat file,
% which can be used by CosInv_v1.0.
%
% The format of input GPS coseismic displacement is
%  Sta   Lon     Lat      E[cm]    N[cm]    U[cm]    SigE     SigN    SigU
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
disp('++++                    Save GPS as Mat File                   ++++');
disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
%% Nepal earthquake
% warning off;clearvars;close all;
% addpath('Nepal Coseismic Disp');
% code_dir='E:\Coseismic\CosInv_v1.0\Code';
% SetPaths(code_dir);
% 
% % Geographic coordinates of epicenter
% origin = [28.147 84.708];
% 
% % [ LOAD COSEISMIC GPS OBSERVATION ]
% disp('    [Load coseismic GPS observations]');
% gps_input_name = 'coseismic_GPS.txt';
% gps_output_name = 'Nepal_GPS';
% gps_file = load(gps_input_name);
% 
% % Number of GPS observation points
% gps.nobs = size(gps_file,1);
% % Geographic position of GPS station
% gps.lon  = gps_file(:,1);
% gps.lat  = gps_file(:,2);
% % Local coordinates of GPS station
% llh_gps  = [gps.lat';gps.lon']; xy_gps  = llh2localxy(llh_gps,origin); 
% gps.x    = xy_gps(:,1);       % East[km]
% gps.y    = xy_gps(:,2);       % North[km]
% % Displacements and variance of GPS station
% gps.disp = gps_file(:,3:5);
% gps.var  = gps_file(:,6:8)'; gps.var = gps.var(:); gps.var = diag(gps.var.^2);
% 
% % [ SAVE THE GPS AS A .MAT FILE ]
% disp('    [Save the GPS as a .MAT file]');
% save(['E:\Coseismic\CosInv_v1.0\Nepal Coseismic Disp\',gps_output_name,'.mat'],'gps');

%% aketao earthquake
% warning off;clearvars;close all;
% addpath('Aketao Coseismic Disp');
% code_dir='E:\Coseismic\CosInv_v1.0\Code';
% SetPaths(code_dir);
% 
% % Geographic coordinates of epicenter
% origin = [39.27 74.04];
% 
% % [ LOAD COSEISMIC GPS OBSERVATION ]
% disp('    [Load coseismic GPS observations]');
% gps_input_name = 'coseismic_GPS.txt';
% gps_output_name = 'Aketao_GPS';
% gps_file = load(gps_input_name);
% 
% % Number of GPS observation points
% gps.nobs = size(gps_file,1);
% % Geographic position of GPS station
% gps.lon  = gps_file(:,1);
% gps.lat  = gps_file(:,2);
% % Local coordinates of GPS station
% llh_gps  = [gps.lat';gps.lon']; xy_gps  = llh2localxy(llh_gps,origin); 
% gps.x    = xy_gps(:,1);       % East[km]
% gps.y    = xy_gps(:,2);       % North[km]
% % Displacements and variance of GPS station
% gps.disp = gps_file(:,3:5);
% gps.var  = gps_file(:,6:8)'; gps.var = gps.var(:); gps.var = diag(gps.var.^2);
% 
% % [ SAVE THE GPS AS A .MAT FILE ]
% disp('    [Save the GPS as a .MAT file]');
% save(['E:\Coseismic\CosInv_v1.0\Aketao Coseismic Disp\',gps_output_name,'.mat'],'gps');

%% Irian earthquake
% warning off;clearvars;close all;
% addpath('Irian Coseismic Disp');
% code_dir='E:\Coseismic\CosInv_v1.0\Code';
% SetPaths(code_dir);
% 
% % Geographic coordinates of epicenter
% origin = [-0.51 132.78];
% 
% % [ LOAD COSEISMIC GPS OBSERVATION ]
% disp('    [Load coseismic GPS observations]');
% gps_input_name = 'Irian_checkerboard_disp_gps.txt';%'coseismic_GPS.txt';
% gps_output_name = 'Irian_checkerboard_GPS';%'Irian_GPS';
% gps_file = load(gps_input_name);
% 
% % Number of GPS observation points
% gps.nobs = size(gps_file,1);
% % Geographic position of GPS station
% gps.lon  = gps_file(:,1);
% gps.lat  = gps_file(:,2);
% % Local coordinates of GPS station
% llh_gps  = [gps.lat';gps.lon']; xy_gps  = llh2localxy(llh_gps,origin); 
% gps.x    = xy_gps(:,1);       % East[km]
% gps.y    = xy_gps(:,2);       % North[km]
% % Displacements and variance of GPS station
% gps.disp = gps_file(:,3:5);
% gps.var  = gps_file(:,6:8)'; gps.var = gps.var(:); gps.var = diag(gps.var.^2);
% 
% % [ SAVE THE GPS AS A .MAT FILE ]
% disp('    [Save the GPS as a .MAT file]');
% save(['E:\Coseismic\CosInv_v1.0\Irian Coseismic Disp\',gps_output_name,'.mat'],'gps');

%% Indonesia earthquake
warning off;clearvars;close all;
addpath('IndonesiaCoseismicDisp');
code_dir='/Users/wangshuai/Downloads/code/CosInv_v1.1/Code';
SetPaths(code_dir);

% Geographic coordinates of epicenter
origin = [-1.328 120.384]; % GCMT 

% [ LOAD COSEISMIC GPS OBSERVATION ]
disp('    [Load coseismic GPS observations]');
gps_input_name = 'coseismic_GPS.txt';%'coseismic_GPS.txt';
gps_output_name = 'Indonesia_GPS';%
gps_file = load(gps_input_name);

% Number of GPS observation points
gps.nobs = size(gps_file,1);
% Geographic position of GPS station
gps.lon  = gps_file(:,1);
gps.lat  = gps_file(:,2);
% Local coordinates of GPS station
llh_gps  = [gps.lat';gps.lon']; xy_gps  = llh2localxy(llh_gps,origin); 
gps.x    = xy_gps(:,1);       % East[km]
gps.y    = xy_gps(:,2);       % North[km]
% Displacements and variance of GPS station
gps.disp = gps_file(:,3:5);
gps.var  = gps_file(:,6:8)'; gps.var = gps.var(:); gps.var = diag(gps.var.^2);

% [ SAVE THE GPS AS A .MAT FILE ]
disp('    [Save the GPS as a .MAT file]');
save(['/Users/wangshuai/Downloads/code/CosInv_v1.1/IndonesiaCoseismicDisp/',gps_output_name,'.mat'],'gps');



