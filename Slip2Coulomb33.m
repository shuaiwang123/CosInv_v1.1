% This script is used to convert the inversion result by CosInv_v1.1 to the
% .inp file required by Coulomb3.3
%
% NOTE
%  1. You must determined which vertices on fault plane are the points of
%    (X_start,Y_start) and (X_fin,Y_fin)
% 
%  2. In COULOMB3.3 right lateral and thrust slip are positive, however, in
%    CosInv_v1.1 left lateral and thrust slip are positive
%
% By Shuai WANG @POLYU 2018-03-29
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
close all;
disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
disp('+++++++++            Welcome To Slip2Coulomb33           ++++++++++');
disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');

%% ------ Poso earthquake 
% % [ 0. Pre-settings]
% disp('    [0. Pre-settings]');
% 
% code_dir='/Users/wangshuai/Downloads/code/CosInv_v1.1/Code';
% SetPaths(code_dir);
% addpath('IndonesiaCoseismicDisp');
% cosinv_filename = 'indonesia_fault_traceMw6.6_Dip44_Southdipping_Len40_Patch1_Smoothing0.03.rec.mat';
% coulomb_inpfile = 'indonesia_cosinv2coulomb.inp';
% origin = [120.384 -1.328];% GCMT 
% 
% % [ 1. Out the heading line ]
% disp('    [1. Out the heading line]');
% 
% fid = fopen(coulomb_inpfile,'w');
% fprintf(fid,'%s\n','CosInv_v1.1 Distributed-SLip Model into Coulomb 3.3');
% fprintf(fid,'%s\n','The fault Parameters from CosInv_v1.1' );
% fprintf(fid,'%s\n',['#reg1=  0  #reg2=  0   #fixed=  ' num2str(n_patches) '  sym=  1']);
% fprintf(fid,'%s\n',' PR1=       .250      PR2=       .250    DEPTH=        0.');
% fprintf(fid,'%s\n','  E1=   0.800000E+06   E2=   0.800000E+06');
% fprintf(fid,'%s\n','XSYM=       .000     YSYM=       .000');
% fprintf(fid,'%s\n','FRIC=       .800');
% fprintf(fid,'%s\n','S1DR=    24.0001     S1DP=      0.0001    S1IN=    100.000     S1GD=   .000000');
% fprintf(fid,'%s\n','S3DR=   114.0001     S3DP=      0.0001    S3IN=     30.000     S3GD=   .000000');
% fprintf(fid,'%s\n','S2DR=    89.9999     S2DP=     -89.999    S2IN=      0.000     S2GD=   .000000');
% 
% 
% % [ 2. Out the fault model ]
% disp('    [2. Out the fault model]');
% 
% fprintf(fid,'%s\n','  #   X-start    Y-start     X-fin      Y-fin   Kode  rt.lat    reverse   dip angle     top      bot');
% fprintf(fid,'%s\n','xxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx');
% 
% for ii = 1:n_patches
%     fprintf(fid,'%3.0f %10.3f %10.3f %10.3f %10.3f %4.0d %10.3f %10.3f %10.2f %10.5f %10.5f %10s\n',...
%         1,   vertices(ii,7),   vertices(ii,8),    vertices(ii,4),   vertices(ii,5),   100,  -slip_x(ii),    slip_y(ii), 44, vertices(ii,9), vertices(ii,3), num2str(ii));
% end
% fprintf(fid,'%s\n','');
% 
% 
% % [ 3. Out grid parameters ]
% disp('    [3. Out grid parameters]');
% 
% % Choose the intresting area from the map
% disp('    ... Choose th interesting area from the map ...');
% figure;
% plot([-100,-100,100,100],[-100,100,100,-100]);hold on;
% scatter(xc,yc,10,'filled');hold on;
% scatter(0,0,100,'rp','filled');grid on;
% 
% hotspot    = ginput;
% start_x    = hotspot(1,1);
% start_y    = hotspot(1,2);
% finish_x   = hotspot(2,1);
% finish_y   = hotspot(2,2);
% 
% hotspot_xy = [hotspot(:,1)';hotspot(:,2)'];
% hotspot_ll = local2llh(hotspot_xy,origin);
% start_lon  = hotspot_ll(1,1);
% start_lat  = hotspot_ll(2,1);
% finish_lon = hotspot_ll(1,2);
% finish_lat = hotspot_ll(2,2);
% 
% 
% fprintf(fid,'%s\n','     Grid Parameters');
% fprintf(fid,'%s\n',['  1  ----------------------------  Start-x =    ' num2str(start_x)]);
% fprintf(fid,'%s\n',['  2  ----------------------------  Start-y =    ' num2str(start_y)]);
% fprintf(fid,'%s\n',['  3  --------------------------   Finish-x =    ' num2str(finish_x)]);
% fprintf(fid,'%s\n',['  4  --------------------------   Finish-y =    ' num2str(finish_y)]);
% fprintf(fid,'%s\n',['  5  ------------------------  x-increment =    ' num2str(2)]);
% fprintf(fid,'%s\n',['  6  ------------------------  y-increment =    ' num2str(2)]);
% 
% 
% % [ 4. Out plot parameters ]
% disp('    [4. Out plot parameters]');
% 
% fprintf(fid,'%s\n','     Size Parameters');
% fprintf(fid,'%s\n','  1  --------------------------  Plot size =        1.0000000');
% fprintf(fid,'%s\n','  2  --------------  Shade/Color increment =        0.2000000');
% fprintf(fid,'%s\n','  3  ------  Exaggeration for disp.& dist. =    10000.0000000');
% 
% 
% % [ 5. Out cross section parameters ]
% disp('    [5. Out cross section parameters]');
% disp('    ... Choose the start/finish points of cross section ...');
% hotspot    = ginput;
% cs_start_x = hotspot(1,1);
% cs_start_y = hotspot(1,2);
% cs_finish_x= hotspot(2,1);
% cs_finish_y= hotspot(2,2);
% 
% fprintf(fid,'%s\n','Cross section default');
% fprintf(fid,'%s\n',['  1  ----------------------------  Start-x =    ' num2str(cs_start_x)]);
% fprintf(fid,'%s\n',['  2  ----------------------------  Start-y =    ' num2str(cs_start_y)]);
% fprintf(fid,'%s\n',['  3  --------------------------   Fiiish-x =    ' num2str(cs_finish_x)]);
% fprintf(fid,'%s\n',['  4  --------------------------   Fiiish-y =    ' num2str(cs_finish_y)]);
% fprintf(fid,'%s\n','  5  ------------------  Distant-increment =     2.000000');
% fprintf(fid,'%s\n','  6  ----------------------------  Z-depth =     30.00000');
% fprintf(fid,'%s\n','  7  ------------------------  Z-increment =     1.000000');
% 
% 
% % [ 6. Out map information ]
% disp('    [6. Out map information]');
% 
% fprintf(fid,'%s\n','     Map information');
% fprintf(fid,'%s\n',['  1  ---------------------------- min. lon =    ' num2str(start_lon)]);
% fprintf(fid,'%s\n',['  2  ---------------------------- max. lon =    ' num2str(finish_lon)]);
% fprintf(fid,'%s\n',['  3  ---------------------------- zero lon =    ' num2str(origin(1))]);
% fprintf(fid,'%s\n',['  4  ---------------------------- min. lat =    ' num2str(start_lat)]);
% fprintf(fid,'%s\n',['  5  ---------------------------- max. lat =    ' num2str(finish_lat)]);
% fprintf(fid,'%s\n',['  6  ---------------------------- zero lat =    ' num2str(origin(2))]);
% fclose(fid);


%% ------ Papua earthquake 
% [ 0. Pre-settings]
disp('    [0. Pre-settings]');

code_dir='/Users/wangshuai/Downloads/code/CosInv_v1.1/Code';
SetPaths(code_dir);
addpath('/Users/wangshuai/Downloads/work/20180226papua/figure/slipJGRrevision');
% cosinv_filename = 'CosInv_papua_FiveSegInTotal_seg34IntoOne_fault_traceMw7.5_dip32.rec.mat';
% cosinv_filename = 'CosInv_papua_seg6_fault_traceMw7.5_dip46_patchsize2km_Lcurve.rec.mat';
cosinv_filename = 'CosInv_papua_postseismic_fault_traceMw7.5_dip25_SlipFreely_Lcurve.rec.mat';
load(cosinv_filename);
% --- Calculate the CFS at depth
% coulomb_inpfile = 'papua_cosinv2coulombSeg12345On6Post.inr';
% coulomb_inpfile = 'tmp2.inp';
%slip_x = slip_x*0;
%slip_y = slip_y*0;
%netslip = Tslip*0;

% netslip           = Tslip;
% netslip(1:255)    = netslip(1:255)*0;
% netslip(811:1170) = netslip(811:1170)*0;



RAKE    = zeros(n_patches,1);
for ii = 1:n_patches
    if slip_x(ii) == 0 & slip_y(ii) ~= 0 
        RAKE(ii) = 90;
    elseif slip_y(ii) == 0 & slip_x(ii) ~= 0
        RAKE(ii) = 0;
    elseif slip_x(ii) ~= 0 & slip_y(ii) ~= 0
        RAKE(ii) = atand(slip_y(ii)/slip_x(ii));
    elseif slip_x(ii) == 0 & slip_y(ii) == 0
        RAKE(ii) = 90;
    end
end

% --- Calculate the CFS on receiver fault 
% coulomb_inpfile = 'papua_cosinv2coulombSection1.inp';
% slip_x(530:end) = 0;
% slip_y(530:end) = 0;

%
% coulomb_inpfile = 'papua_cosinv2coulombSection12.inp';
% slip_x(714:end) = 0;
% slip_y(714:end) = 0;

%
% coulomb_inpfile = 'papua_cosinv2coulombSection123.inp';
% slip_x(921:end) = 0;
% slip_y(921:end) = 0;

%
% coulomb_inpfile = 'papua_cosinv2coulombSection1234.inp';


origin = [142.754 -6.070 ]; % USGS

% [ 1. Out the heading line ]
disp('    [1. Out the heading line]');

fid = fopen(coulomb_inpfile,'w');
fprintf(fid,'%s\n','CosInv_v1.1 Distributed-SLip Model into Coulomb 3.3');
fprintf(fid,'%s\n','The fault Parameters from CosInv_v1.1' );
fprintf(fid,'%s\n',['#reg1=  0  #reg2=  0   #fixed=  ' num2str(n_patches) '  sym=  1']);
fprintf(fid,'%s\n',' PR1=       .250      PR2=       .250    DEPTH=        0.');
fprintf(fid,'%s\n','  E1=   0.800000E+06   E2=   0.800000E+06');
fprintf(fid,'%s\n','XSYM=       .000     YSYM=       .000');
fprintf(fid,'%s\n','FRIC=       .800');
fprintf(fid,'%s\n','S1DR=    24.0001     S1DP=      0.0001    S1IN=    100.000     S1GD=   .000000');
fprintf(fid,'%s\n','S3DR=   114.0001     S3DP=      0.0001    S3IN=     30.000     S3GD=   .000000');
fprintf(fid,'%s\n','S2DR=    89.9999     S2DP=     -89.999    S2IN=      0.000     S2GD=   .000000');


% [ 2. Out the fault model ]
disp('    [2. Out the fault model]');

% --- Out in strike- and dip- components format
% fprintf(fid,'%s\n','  #   X-start    Y-start     X-fin      Y-fin   Kode  rt.lat    reverse   dip angle     top      bot');
% fprintf(fid,'%s\n','xxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx');
% 
% for ii = 1:n_patches
%     fprintf(fid,'%3.0f %10.3f %10.3f %10.3f %10.3f %4.0d %10.3f %10.3f %10.2f %10.5f %10.5f %10s\n',...
%         1,   vertices(ii,7),   vertices(ii,8),    vertices(ii,4),   vertices(ii,5),   100,  -slip_x(ii),    slip_y(ii), dip(ii), vertices(ii,9), vertices(ii,3), num2str(ii));
% end
% fprintf(fid,'%s\n','');

% --- Out in netslip format
fprintf(fid,'%s\n','  #   X-start    Y-start     X-fin      Y-fin   Kode  rake    netslip   dip angle     top      bot');
fprintf(fid,'%s\n','xxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx');

for ii = 1:n_patches
    fprintf(fid,'%3.0f %10.3f %10.3f %10.3f %10.3f %4.0d %10.3f %10.3f %10.2f %10.5f %10.5f %10s\n',...
        1,   vertices(ii,7),   vertices(ii,8),    vertices(ii,4),   vertices(ii,5),   100,  RAKE(ii),    netslip(ii), dip(ii), vertices(ii,9), vertices(ii,3), num2str(ii));
end
fprintf(fid,'%s\n','');

% [ 3. Out grid parameters ]
disp('    [3. Out grid parameters]');

% Choose the intresting area from the map
disp('    ... Choose th interesting area from the map ...');
figure;
plot([-100,-100,100,100],[-100,100,100,-100]);hold on;
scatter(xc,yc,10,'filled');hold on;
scatter(0,0,100,'rp','filled');grid on;

hotspot    = ginput;
start_x    = hotspot(1,1);
start_y    = hotspot(1,2);
finish_x   = hotspot(2,1);
finish_y   = hotspot(2,2);

hotspot_xy = [hotspot(:,1)';hotspot(:,2)'];
hotspot_ll = local2llh(hotspot_xy,origin);
start_lon  = hotspot_ll(1,1);
start_lat  = hotspot_ll(2,1);
finish_lon = hotspot_ll(1,2);
finish_lat = hotspot_ll(2,2);


fprintf(fid,'%s\n','     Grid Parameters');
fprintf(fid,'%s\n',['  1  ----------------------------  Start-x =    ' num2str(start_x)]);
fprintf(fid,'%s\n',['  2  ----------------------------  Start-y =    ' num2str(start_y)]);
fprintf(fid,'%s\n',['  3  --------------------------   Finish-x =    ' num2str(finish_x)]);
fprintf(fid,'%s\n',['  4  --------------------------   Finish-y =    ' num2str(finish_y)]);
fprintf(fid,'%s\n',['  5  ------------------------  x-increment =    ' num2str(2)]);
fprintf(fid,'%s\n',['  6  ------------------------  y-increment =    ' num2str(2)]);


% [ 4. Out plot parameters ]
disp('    [4. Out plot parameters]');

fprintf(fid,'%s\n','     Size Parameters');
fprintf(fid,'%s\n','  1  --------------------------  Plot size =        1.0000000');
fprintf(fid,'%s\n','  2  --------------  Shade/Color increment =        0.2000000');
fprintf(fid,'%s\n','  3  ------  Exaggeration for disp.& dist. =    10000.0000000');


% [ 5. Out cross section parameters ]
disp('    [5. Out cross section parameters]');
disp('    ... Choose the start/finish points of cross section ...');
hotspot    = ginput;
cs_start_x = hotspot(1,1);
cs_start_y = hotspot(1,2);
cs_finish_x= hotspot(2,1);
cs_finish_y= hotspot(2,2);

fprintf(fid,'%s\n','Cross section default');
fprintf(fid,'%s\n',['  1  ----------------------------  Start-x =    ' num2str(cs_start_x)]);
fprintf(fid,'%s\n',['  2  ----------------------------  Start-y =    ' num2str(cs_start_y)]);
fprintf(fid,'%s\n',['  3  --------------------------   Fiiish-x =    ' num2str(cs_finish_x)]);
fprintf(fid,'%s\n',['  4  --------------------------   Fiiish-y =    ' num2str(cs_finish_y)]);
fprintf(fid,'%s\n','  5  ------------------  Distant-increment =     2.000000');
fprintf(fid,'%s\n','  6  ----------------------------  Z-depth =     60.00000');
fprintf(fid,'%s\n','  7  ------------------------  Z-increment =     1.000000');


% [ 6. Out map information ]
disp('    [6. Out map information]');

fprintf(fid,'%s\n','     Map information');
fprintf(fid,'%s\n',['  1  ---------------------------- min. lon =    ' num2str(start_lon)]);
fprintf(fid,'%s\n',['  2  ---------------------------- max. lon =    ' num2str(finish_lon)]);
fprintf(fid,'%s\n',['  3  ---------------------------- zero lon =    ' num2str(origin(1))]);
fprintf(fid,'%s\n',['  4  ---------------------------- min. lat =    ' num2str(start_lat)]);
fprintf(fid,'%s\n',['  5  ---------------------------- max. lat =    ' num2str(finish_lat)]);
fprintf(fid,'%s\n',['  6  ---------------------------- zero lat =    ' num2str(origin(2))]);
fclose(fid);