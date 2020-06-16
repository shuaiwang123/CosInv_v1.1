% Script to generate the noise datasets using the Monte Carlo analysis,which
% can be used to estimate the slip distribution uncertainties 
%
% REF.
% B. Parsonset et al, 2006, GJI. The 1994 Sefidabeh (eastern Iran) earthquakes 
%    revisited: new evidence from satellite radar interferometry and carbonate 
%    dating about the growth of an active fold above a blind thrust fault. 
%    
% By shwang @whu 2017-03-12
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clear;close all;clearvars;

%% ==== Aketao earthquake ====
% 
% addpath('Aketao Coseismic Disp');
% origin   = [74.09 39.27];
% code_dir ='E:\Coseismic\CosInv_v1.0\Code';
% SetPaths(code_dir);
% 
% %% [Load InSAR data ]
% insar_input_name = 'Aketao_InSAR2.mat';
% 
% load(insar_input_name);
% % number of InSAR pixels
% nobs_mc = insar.nobs;      
% lon_mc  = insar.lon;
% lat_mc  = insar.lat;
% x_mc    = insar.x;
% y_mc    = insar.y;
% vect_mc = insar.vect;
% los_mc  = insar.los;
% var_mc  = insar.var;
% 
% L       = chol(var_mc);
% 
% %% [ Random number created by Monte Carlo ]
% % the average and standard deviation used for Monte Carlo
% mu    = 0;
% sigma = 1;
% % number of dataset to be created using Monte Carlo
% num_mc_datasets = 100;
% 
% for ii=1:num_mc_datasets
%     % Guassian uncorrelated noise, with a mean of mu and a standard
%     % deviation sigma, please refer to Parsons 2006 GJI
%     xx_mc = normrnd(mu,sigma,nobs_mc,1);
%     yy_mc = L*xx_mc; 
%     los_mc    = los_mc + yy_mc;
%     
%     insar.nobs = nobs_mc;
%     insar.lon  = lon_mc;
%     insar.lat  = lat_mc;
%     insar.x    = x_mc;
%     insar.y    = y_mc;
%     insar.vect = vect_mc;
%     insar.los  = los_mc;
%     insar.var  = var_mc;
%     
%     % save result
%     insar_output_name = ['mc_',num2str(ii),'_',insar_input_name];
%     save(['E:\Coseismic\CosInv_v1.0\Aketao Coseismic Disp\montecarlo\',insar_output_name],'insar');
% end

% %% ==== Irian earthquake ====
% 
% addpath('Irian Coseismic Disp');
% origin = [132.78 -0.51];
% code_dir ='E:\Coseismic\CosInv_v1.0\Code';
% SetPaths(code_dir);
% 
% %% [Load InSAR data ]
% insar_input_name = 'Irian_InSAR3.mat';
% 
% load(insar_input_name);
% % number of InSAR pixels
% nobs_mc = insar.nobs;      
% lon_mc  = insar.lon;
% lat_mc  = insar.lat;
% x_mc    = insar.x;
% y_mc    = insar.y;
% vect_mc = insar.vect;
% los_mc  = insar.los;
% var_mc  = insar.var;
% 
% L       = chol(var_mc);
% 
% %% [ Random number created by Monte Carlo ]
% % the average and standard deviation used for Monte Carlo
% mu    = 0;
% sigma = 1;
% % number of dataset to be created using Monte Carlo
% num_mc_datasets = 100;
% 
% for ii=1:num_mc_datasets
%     % Guassian uncorrelated noise, with a mean of mu and a standard
%     % deviation sigma, please refer to Parsons 2006 GJI
%     xx_mc = normrnd(mu,sigma,nobs_mc,1);
%     yy_mc = L*xx_mc; 
%     los_mc    = los_mc + yy_mc;
%     
%     insar.nobs = nobs_mc;
%     insar.lon  = lon_mc;
%     insar.lat  = lat_mc;
%     insar.x    = x_mc;
%     insar.y    = y_mc;
%     insar.vect = vect_mc;
%     insar.los  = los_mc;
%     insar.var  = var_mc;
%     
%     % save result
%     insar_output_name = ['mc_',num2str(ii),'_',insar_input_name];
%     save(['E:\Coseismic\CosInv_v1.0\Irian Coseismic Disp\montecarlo\',insar_output_name],'insar');
% end


%% ==== Indonesia earthquake ====

addpath('IndonesiaCoseismicDisp');
origin = [120.384 -1.328];% GCMT 
code_dir ='/Users/wangshuai/Downloads/code/CosInv_v1.1';
SetPaths(code_dir);

%% [Load InSAR data ]
insar_input_name = 'Indonesia_InSAR1.mat';

load(insar_input_name);
% number of InSAR pixels
nobs_mc = insar.nobs;      
lon_mc  = insar.lon;
lat_mc  = insar.lat;
x_mc    = insar.x;
y_mc    = insar.y;
vect_mc = insar.vect;
LOS_mc  = insar.los;
var_mc  = insar.var;

L       = chol(var_mc);

%% [ Random number created by Monte Carlo ]
% the average and standard deviation used for Monte Carlo
mu    = 0;
% --- dataset 1
sigma = 1;%2.33;
% --- dataset 2
% sigma = 1;%1.7734;
% number of dataset to be created using Monte Carlo
num_mc_datasets = 100;

for ii=1:num_mc_datasets
    % Guassian uncorrelated noise, with a mean of mu and a standard
    % deviation sigma, please refer to Parsons 2006 GJI
    xx_mc = normrnd(mu,sigma,nobs_mc,1);
    %index = find(abs(xx_mc)>=3);
    %xx_mc(index) = xx_mc(index)/max(abs(xx_mc));
    yy_mc = xx_mc;%L*xx_mc; 
    los_mc= LOS_mc + yy_mc;
    
    insar.nobs = nobs_mc;
    insar.lon  = lon_mc;
    insar.lat  = lat_mc;
    insar.x    = x_mc;
    insar.y    = y_mc;
    insar.vect = vect_mc;
    insar.los  = los_mc;
    insar.var  = var_mc;
    
    % save result
    insar_output_name = ['mc_',num2str(ii),'_',insar_input_name];
    save(['/Users/wangshuai/Downloads/code/CosInv_v1.1/IndonesiaCoseismicDisp/montecarlo/',insar_output_name],'insar');
end

% %% ==== Australia earthquake ====
% 
% addpath('AustraliaCoseismicDisp');
% origin = [129.884 -25.566];% GCMT 
% code_dir ='/Users/wangshuai/Downloads/code/CosInv_v1.1';
% SetPaths(code_dir);
% 
% %% [Load InSAR data ]
% insar_input_name = 'Aus_InSAR2.mat';
% 
% load(insar_input_name);
% % number of InSAR pixels
% nobs_mc = insar.nobs;      
% lon_mc  = insar.lon;
% lat_mc  = insar.lat;
% x_mc    = insar.x;
% y_mc    = insar.y;
% vect_mc = insar.vect;
% LOS_mc  = insar.los;
% var_mc  = insar.var;
% 
% L       = chol(var_mc);
% 
% %% [ Random number created by Monte Carlo ]
% % the average and standard deviation used for Monte Carlo
% mu    = 0;
% % --- dataset 1
% % sigma = 0.9110;
% % --- dataset 2
% sigma = 0.9165;
% % number of dataset to be created using Monte Carlo
% num_mc_datasets = 100;
% 
% for ii=1:num_mc_datasets
%     % Guassian uncorrelated noise, with a mean of mu and a standard
%     % deviation sigma, please refer to Parsons 2006 GJI
%     xx_mc = normrnd(mu,sigma,nobs_mc,1);
%     %index = find(abs(xx_mc)>=3);
%     %xx_mc(index) = xx_mc(index)/max(abs(xx_mc));
%     yy_mc = xx_mc;%L*xx_mc; 
%     los_mc= LOS_mc + yy_mc;
%     
%     insar.nobs = nobs_mc;
%     insar.lon  = lon_mc;
%     insar.lat  = lat_mc;
%     insar.x    = x_mc;
%     insar.y    = y_mc;
%     insar.vect = vect_mc;
%     insar.los  = los_mc;
%     insar.var  = var_mc;
%     
%     % save result
%     insar_output_name = ['mc_',num2str(ii),'_',insar_input_name];
%     save(['/Users/wangshuai/Downloads/code/CosInv_v1.1/AustraliaCoseismicDisp/montecarlo/',insar_output_name],'insar');
% end