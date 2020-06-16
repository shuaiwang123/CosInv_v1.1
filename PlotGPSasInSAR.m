% PlotGPSasInSAR.M
%
% Script to project the GPS measurements into LOS direction and plot them.
%
% by shwang @whu 2017-05-22
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%==========================================================================
%======                                                             =======
%===========                     Pre Setting               ================
%==========================================================================
clc;clearvars;close all;
disp('...... Pre Setting ......');
addpath('Irian Coseismic Disp');
origin = [132.78 -0.51];
code_dir='E:\Coseismic\CosInv_v1.0\Code';
SetPaths(code_dir);

% Convert units
m2cm           = 100;

% Epicenter
epicenter      = [132.78,-0.51;133.34,-0.71];

% GPS and InSAR file
gps_file_name   = 'E:\Work\Irian\p000gs2d.disp.txt';
% insar_file_name = 'E:\Coseismic\CosInv_v1.0\Irian Coseismic Disp\ALOS_T386A_Gamma.QT';
% insar_file_name = 'E:\Coseismic\CosInv_v1.0\Irian Coseismic Disp\ALOS_T387A_Gamma.QT';
 insar_file_name = 'E:\Coseismic\CosInv_v1.0\Irian Coseismic Disp\ALOS_T057D_Gamma.QT';

% Mask file
mask_file_name  = 'E:\Work\Irian\irian_coast_shore.dat';

%==========================================================================
%======                                                             =======
%========  Load GPS and InSAR[including looking vector] Files    ==========
%==========================================================================
disp('...... Load GPS and InSAR Files ......');
disp('    [Load GPS File]');
if ~exist(gps_file_name,'file'),
    error('    [Error: GPS file does not exist,plese check!!!]');
else
    gps_file = load(gps_file_name);
    gps_lon  = gps_file(:,1);
    gps_lat  = gps_file(:,2);
    gps_e    = gps_file(:,3)*m2cm;
    gps_n    = gps_file(:,4)*m2cm;
    gps_u    = gps_file(:,5)*m2cm;
end

disp('    [Load InSAR File]');
if ~exist(insar_file_name,'file'),
    error('    [Error: InSAR file does not exist,plese check!!!]');
else
    insar_file = load(insar_file_name);
    insar_lon  = insar_file(:,1);
    insar_lat  = insar_file(:,2);
    insar_lvE  = insar_file(:,3);
    insar_lvN  = insar_file(:,4);
    insar_lvU  = insar_file(:,5);
    insar_los  = insar_file(:,6);
end

%==========================================================================
%======                                                             =======
%========              Project GPS onto LOS Direction            ==========
%==========================================================================
disp('...... Project GPS into LOS Direction ......');

disp('    [Interpolate Looking Vector on GPS Sites]');
gps_lvE = griddata(insar_lon,insar_lat,insar_lvE,gps_lon,gps_lat,'v4');
gps_lvN = griddata(insar_lon,insar_lat,insar_lvN,gps_lon,gps_lat,'v4');
gps_lvU = griddata(insar_lon,insar_lat,insar_lvU,gps_lon,gps_lat,'v4');

disp('    [Project GPS deformation onto LOS Direction]');
gps_los = gps_e.*gps_lvE + gps_n.*gps_lvN + gps_u.*gps_lvU;


%==========================================================================
%======                                                             =======
%========        Extract Stations and Pixels by Mask File        ==========
%==========================================================================
% disp('...... Extract Stations and Pixels by Mask File ......');
% 
% mask_file = load(mask_file_name);
% mask_lon  = mask_file(:,1);
% mask_lat  = mask_file(:,2);
% 
% disp('    [Extract Stations]');
% gps_in    = inpolygon(gps_lon,gps_lat,mask_lon,mask_lat);
% gps_index = find(gps_in ~= 0);
% 
% gps_lon   = gps_lon(gps_index);
% gps_lat   = gps_lat(gps_index);
% gps_e     = gps_e(gps_index);
% gps_n     = gps_n(gps_index);
% gps_u     = gps_u(gps_index);
% gps_los   = gps_los(gps_index);


%==========================================================================
%======                                                             =======
%========                         Plotting                       ==========
%==========================================================================
disp('...... Plotting ......');

figure;
% [ Plot GPS-type Measurements]
disp('    [Plot GPS-type Measurements]');
subplot(1,3,1);
% Plot vertical deformation
plotShadeMap(gps_lon,gps_lat,gps_u,mask_file_name,'GPS-type Deformation');hold on;
% Plot horizontal deformation
quiver(gps_lon,gps_lat,gps_e,gps_n,'w');
% Plot epicenter
for ii = 1:size(epicenter,1);
    scatter(epicenter(ii,1),epicenter(ii,2),100,'p','w','filled');
end
% xlim,ylim
xlim([131,134]);
ylim([-1.5,0]);
axis square;

% [ Plot GPS Projected Deformation ]
disp('    [Plot the GPS Projected Deformation[LOS]');
subplot(1,3,2);
hold on;
plotShadeMap(gps_lon,gps_lat,gps_los,mask_file_name,'GPS-type Projected Deformation');
% Plot epicenter
for ii = 1:size(epicenter,1);
    scatter(epicenter(ii,1),epicenter(ii,2),100,'p','w','filled');
end
% xlim,ylim
xlim([131,134]);
ylim([-1.5,0]);
axis square;

% [ Plot the INSAR Measurements ]
disp('    [Plot the InSAR Measurements');
subplot(1,3,3);
hold on;
plotShadeMap(insar_lon,insar_lat,insar_los,mask_file_name,'InSAR Measurements');
% Plot epicenter
for ii = 1:size(epicenter,1);
    scatter(epicenter(ii,1),epicenter(ii,2),100,'p','w','filled');
end
% xlim,ylim
xlim([132,134]);
ylim([-1.5,0]);
axis square;

