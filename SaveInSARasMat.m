% This script is used to save the InSAR coseismic LOS displacement as .mat
% file, which can be used by CosInv_v1.0.
%
% The format of input InSAR coseismic LOS displacement is
%   LON       LAT      E_vector     N_vector     U_vector    LOS[cm]
% 
% ### Note that no matter what the original unit of the observed LOS data,
%     we should convert it into the unit required by CosInv_v1.0.
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
disp('++++                   Save InSAR as Mat File                  ++++');
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
% % [ LOAD COSEISMIC InSAR OBSERVATION ]
% disp('    [Load coseismic InSAR observations]');
% % ... InSAR dataset 1 ...
% insar_input_name = 'LOS2_P48_20150222_20150503_ScanSAR.QT';
% insar_output_name = 'Nepal_InSAR1';
% % LOS variance[cm^2];
% var = 1.44;
% a   = 0.43;
% b   = 10.74;
% % ... InSAR dataset 2 ...
% % insar_input_name = 'LOS2_P157_20150221_20150502_StriMap.QT';
% % insar_output_name = 'Nepal_InSAR2';
% % % LOS variance[cm^2];
% % var = 1.23;
% % a   = 0.47;
% % b   = 16.25;
% 
% insar_file = load(insar_input_name);
% 
% % Number of InSAR pixels
% insar.nobs = size(insar_file,1);
% % Geographic position of InSAR pixel
% insar.lon  = insar_file(:,1);
% insar.lat  = insar_file(:,2);
% % Local coordinate of InSAR pixel
% llh_insar  = [insar.lat';insar.lon']; xy_insar = llh2localxy(llh_insar,origin);
% insar.x    = xy_insar(:,1);   % East[km]
% insar.y    = xy_insar(:,2);   % North[km]
% % LOS and variance of InSAR pixel
% insar.vect = insar_file(:,3:5);
% insar.los  = insar_file(:,6);
% 
% % [ CALCULATE CONVARIANCE BETWEEN PIXELS ]
% disp('    [Calculate convariance between pixels using empirical fuction]');
% insar.var  = zeros(insar.nobs,insar.nobs);
% distance   = zeros(insar.nobs,insar.nobs);
% for ii=1:insar.nobs
%     for kk=ii:insar.nobs
%         distance(ii,kk) = dist(insar.x(kk),insar.y(kk),insar.x(ii),insar.y(ii));
%     end
% end
% distance = distance + distance';
% 
% for jj=1:insar.nobs
%     % convariance between pixels
%     insar.var(jj,:)  = a*exp(-distance(jj,:)./b);
%     % variance 
%     insar.var(jj,jj) = var;
% end
% 
% % [ SAVE THE INSAR AS A .MAT FILE]
% disp('    [Save InSAR as a .MAT file]');
% save(['E:\Coseismic\CosInv_v1.0\Nepal Coseismic Disp\',insar_output_name,'.mat'],'insar');

%% Aketao earthquake
% warning off;clearvars;close all;
% addpath('Aketao Coseismic Disp');
% code_dir='E:\Coseismic\CosInv_v1.0\Code';
% SetPaths(code_dir);
% 
% % Geographic coordinates of epicenter
% origin = [39.27 74.04];
% 
% % [ LOAD COSEISMIC InSAR OBSERVATION ]
% disp('    [Load coseismic InSAR observations]');
% % % ... InSAR dataset 1 ...
% % insar_input_name = 'S1A_P027A_20161113_20161207_Gamma.QT';
% % insar_output_name = 'Aketao_InSAR1';
% % % LOS variance[cm^2];
% % var = 0.67;
% % a   = 0.67;
% % b   = 5.3;
% % ... InSAR dataset 2 ...
% % insar_input_name = 'S1B_P107D_20161101_20161219_Gamma.QT';
% % insar_output_name = 'Aketao_InSAR2';
% % % LOS variance[cm^2];
% % var = 0.98;
% % a   = 0.98;
% % b   = 9;
% ... InSAR dataset 3 ...
% insar_input_name = 'ALOS2_P162A_20160720_20161207_Gamma.QT';
% insar_output_name = 'Aketao_InSAR3';
% % LOS variance[cm^2];
% var = 0.49;
% a   = 0.49;
% b   = 10;
% 
% insar_file = load(insar_input_name);
% 
% % Number of InSAR pixels
% insar.nobs = size(insar_file,1);
% % Geographic position of InSAR pixel
% insar.lon  = insar_file(:,1);
% insar.lat  = insar_file(:,2);
% % Local coordinate of InSAR pixel
% llh_insar  = [insar.lat';insar.lon']; xy_insar = llh2localxy(llh_insar,origin);
% insar.x    = xy_insar(:,1);   % East[km]
% insar.y    = xy_insar(:,2);   % North[km]
% % LOS and variance of InSAR pixel
% insar.vect = insar_file(:,3:5);
% insar.los  = insar_file(:,6);
% 
% % [ CALCULATE CONVARIANCE BETWEEN PIXELS ]
% disp('    [Calculate convariance between pixels using empirical fuction]');
% insar.var  = zeros(insar.nobs,insar.nobs);
% distance   = zeros(insar.nobs,insar.nobs);
% for ii=1:insar.nobs
%     for kk=ii:insar.nobs
%         distance(ii,kk) = dist(insar.x(kk),insar.y(kk),insar.x(ii),insar.y(ii));
%     end
% end
% distance = distance + distance';
% 
% for jj=1:insar.nobs
%     % convariance between pixels
%     insar.var(jj,:)  = a*exp(-distance(jj,:)./b);
%     % variance 
%     insar.var(jj,jj) = var;
% end
% 
% % [ SAVE THE INSAR AS A .MAT FILE]
% disp('    [Save InSAR as a .MAT file]');
% save(['E:\Coseismic\CosInv_v1.0\Aketao Coseismic Disp\',insar_output_name,'.mat'],'insar');

%% Irian earthquake
% warning off;clearvars;close all;
% addpath('Irian Coseismic Disp');
% code_dir='E:\Coseismic\CosInv_v1.0\Code';
% SetPaths(code_dir);
% 
% % Geographic coordinates of epicenter
% origin = [-0.51 132.78];
% 
% % [ LOAD COSEISMIC InSAR OBSERVATION ]
% disp('    [Load coseismic InSAR observations]');
% % % ... InSAR dataset 1 ...
% % insar_input_name = 'ALOS_T057D_Gamma.QT';
% % insar_output_name = 'Irian_InSAR1';
% % % LOS variance[cm^2];
% % var = 10.11;
% % a   = 10.11;
% % b   = 5.2;
% % ... InSAR dataset 2 ...
% % insar_input_name = 'ALOS_T386A_Gamma.QT';
% % insar_output_name = 'Irian_InSAR2';
% 
% % insar_input_name = 'ALOS_T386A_Gamma_Mw7.6CosEffRev.QT';
% % insar_output_name = 'Irian_InSARCosEffRev2';
% 
% % insar_input_name = 'ALOS_T386A_Gamma_Mw7.6CosEffRev_LinearRampRev.QT';
% % insar_output_name = 'Irian_InSARCosEffRev2_LinearRampRev';
% 
% insar_input_name = 'Irian_checkerboard_disp_insar1.txt';%'ALOS_T386A_Gamma_LinearRampRev.QT';
% insar_output_name = 'Irian_checkerboard_t386';%'Irian_InSARLinearRampRev2';
% % LOS variance[cm^2];
% var = 2.64;
% a   = 2.64;
% b   = 21.82;
% 
% % ... InSAR dataset 3 ...
% % insar_input_name = 'Irian_checkerboard_disp_insar2.txt';%'ALOS_T387A_Gamma.QT';
% % insar_output_name = 'Irian_checkerboard_t387';%'Irian_InSAR3';
% % % LOS variance[cm^2];
% % var = 3.83;
% % a   = 3.83;
% % b   = 4;
% 
% insar_file = load(insar_input_name);
% 
% % Number of InSAR pixels
% insar.nobs = size(insar_file,1);
% % Geographic position of InSAR pixel
% insar.lon  = insar_file(:,1);
% insar.lat  = insar_file(:,2);
% % Local coordinate of InSAR pixel
% llh_insar  = [insar.lat';insar.lon']; xy_insar = llh2localxy(llh_insar,origin);
% insar.x    = xy_insar(:,1);   % East[km]
% insar.y    = xy_insar(:,2);   % North[km]
% % LOS and variance of InSAR pixel
% insar.vect = insar_file(:,3:5);
% insar.los  = insar_file(:,6);
% 
% % [ CALCULATE CONVARIANCE BETWEEN PIXELS ]
% disp('    [Calculate convariance between pixels using empirical fuction]');
% insar.var  = zeros(insar.nobs,insar.nobs);
% distance   = zeros(insar.nobs,insar.nobs);
% for ii=1:insar.nobs
%     for kk=ii:insar.nobs
%         distance(ii,kk) = dist(insar.x(kk),insar.y(kk),insar.x(ii),insar.y(ii));
%     end
% end
% distance = distance + distance';
% 
% for jj=1:insar.nobs
%     % convariance between pixels
%     insar.var(jj,:)  = a*exp(-distance(jj,:)./b);
%     % variance 
%     insar.var(jj,jj) = var;
% end
% 
% % [ SAVE THE INSAR AS A .MAT FILE]
% disp('    [Save InSAR as a .MAT file]');
% save(['E:\Coseismic\CosInv_v1.0\Irian Coseismic Disp\',insar_output_name,'.mat'],'insar');

%% Indonesia earthquake
warning off;clearvars;close all;
addpath('IndonesiaCoseismicDisp');
code_dir='/Users/wangshuai/Downloads/code/CosInv_v1.1/Code';
SetPaths(code_dir);

% Geographic coordinates of epicenter
origin = [-1.328 120.384]; % GCMT 

% [ LOAD COSEISMIC InSAR OBSERVATION ]
disp('    [Load coseismic InSAR observations]');
% % ... InSAR dataset 1 ...
% insar_input_name = 'S1A_P027A_20161113_20161207_Gamma.QT';
% insar_output_name = 'Aketao_InSAR1';
% % LOS variance[cm^2];
% var = 0.67;
% a   = 0.67;
% b   = 5.3;
% ... InSAR dataset 2 ...
% insar_input_name = 'S1B_P107D_20161101_20161219_Gamma.QT';
% insar_output_name = 'Aketao_InSAR2';
% % LOS variance[cm^2];
% var = 0.98;
% a   = 0.98;
% b   = 9;
... InSAR dataset 3 ...
insar_input_name = 'ALOS2_P162A_20160720_20161207_Gamma.QT';
insar_output_name = 'Aketao_InSAR3';
% LOS variance[cm^2];
var = 0.49;
a   = 0.49;
b   = 10;

insar_file = load(insar_input_name);

% Number of InSAR pixels
insar.nobs = size(insar_file,1);
% Geographic position of InSAR pixel
insar.lon  = insar_file(:,1);
insar.lat  = insar_file(:,2);
% Local coordinate of InSAR pixel
llh_insar  = [insar.lat';insar.lon']; xy_insar = llh2localxy(llh_insar,origin);
insar.x    = xy_insar(:,1);   % East[km]
insar.y    = xy_insar(:,2);   % North[km]
% LOS and variance of InSAR pixel
insar.vect = insar_file(:,3:5);
insar.los  = insar_file(:,6);

% [ CALCULATE CONVARIANCE BETWEEN PIXELS ]
disp('    [Calculate convariance between pixels using empirical fuction]');
insar.var  = zeros(insar.nobs,insar.nobs);
distance   = zeros(insar.nobs,insar.nobs);
for ii=1:insar.nobs
    for kk=ii:insar.nobs
        distance(ii,kk) = dist(insar.x(kk),insar.y(kk),insar.x(ii),insar.y(ii));
    end
end
distance = distance + distance';

for jj=1:insar.nobs
    % convariance between pixels
    insar.var(jj,:)  = a*exp(-distance(jj,:)./b);
    % variance 
    insar.var(jj,jj) = var;
end

% [ SAVE THE INSAR AS A .MAT FILE]
disp('    [Save InSAR as a .MAT file]');
save(['E:\Coseismic\CosInv_v1.0\Aketao Coseismic Disp\',insar_output_name,'.mat'],'insar');











