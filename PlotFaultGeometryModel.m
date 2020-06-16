% Usage: script to give a snapshot of fault model
%
% By Shuai WANG @WHU 2018-08-19
%==========================================================================
%clear;close all;
% addpath('PapuaCoseismicDisp');

% [ 0. Pre-settings ]
% ------ Papua earthquake
% fault_model_file = 'papua_FiveSegInTotal_fault_traceMw7.5_dip18.rec';
% fault_model_file = 'papua_FiveSegInTotal_seg34IntoOne_fault_traceMw7.5_dip32.rec';
% fault_model_file = 'papua_RampFlat_fault_traceMw7.5_dip40-18_Depth9.rec';
% fault_model_file = 'papua_SixSegInTotal_fault_traceMw7.5_dip18.rec';
% fault_model_file = 'papua_seg6_fault_traceMw7.5_dip46_patchsize2km.rec';
% fault_model_file = 'papua_postseismic_fault_traceMw7.5_dip25.rec';

% ------ Aketao earthquake
% fault_model_file = 'AketaoCoseismicDisp/GRL/muji_faultTraceSeg123_78.rec';
fault_model_file = 'ChangningCoseisDisp/changning_faultTrace_aftershocks_northeastDipping10.rec';

fault_model      = load(fault_model_file);
%fault_model = fault_model(1110:end,:);

% [ 1. Options for Plot ]
% Define the view angles using median strike and dip
AZ      = 10;
EL      = 40;
ViewAngles = [AZ,EL];

% Defining the color map
% mycolormap = colormap(flipud(hot));
mycolormap = colormap(gray);
% mycolormap = colormap(flipud(gray));

% SLip and slip vector
n_patches = size(fault_model,1);
Tslip = zeros(n_patches,1);

slip_vector = repmat([1 0 0],n_patches,1);

% Options to plot the slip vector
vector_scale = 0.0;
vector_color = 'w';
vector_width = 1;
% options to plot GPS and patch number 
gps_patchnumber = 1;%input('    [Please choosing whether plotting the station name and disp vectors of GPS and the patch number: 0. NO  1.YES]');

% Various options for plot_rectangle_slip_distribution_and_station_vector 
plot_options = {...
    'Perspective',ViewAngles,...
    'FieldVector',slip_vector,vector_scale,vector_color,vector_width,...
    'ColorMap',mycolormap,...
    'AutoScale',1,...
    'Gps_patchnumber',gps_patchnumber,...
    'Figure',...
    'Old',...
    'ColorBar',...
    'Shading','Flat',...
    'ColorBarLabel',['Slip (m)']};

plot_rec_slip_dis(fault_model,Tslip,plot_options);

