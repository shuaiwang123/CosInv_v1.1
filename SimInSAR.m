% This scripts is used to simulate the LOS displacement/unwrapping phase/wrapping phase 
% based on the inverted slip distribution model.
% 
% Format of input InSAR coseismic LOS displacement
%   LON       LAT      E_vector     N_vector     U_vector    LOS[m]
%
% 2017-02-04 By ShuaiWang @Zhengzhou
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
disp('++++++++++++++++        Welcome To SimInSAR        ++++++++++++++++');
disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
addpath('Aketao Coseismic Disp');
origin   = [74.09 39.27];
code_dir ='E:\Coseismic\CosInv_v1.0\Code';
SetPaths(code_dir);

%% 1.Load InSAR data 
disp(['...... [Load InSAR data] ......']);
% === Aketao earthquake ===
% --- T027A ---
insar_obs_file  = 'S1A_P027A_20161113_20161207_Gamma.txt';
insar_sim_file  = 'S1A_P027A_20161113_20161207_Gamma_sim.txt';
% C band(m) 0.0555041577
% L band(m) 0.236
wavelength      = 0.0555041577; 
% % --- T107D --- 
% insar_obs_file  = 'S1B_P107D_20161101_20161219_Gamma.txt';
% insar_sim_file  = 'S1B_P107D_20161101_20161219_Gamma_sim.txt';
% % C band(m) 0.0555041577
% % L band(m) 0.236
% wavelength      = 0.0555041577; 
% % --- T162A --- 
% insar_obs_file  = 'ALOS2_P162A_20160720_20161207_Gamma.txt';
% insar_sim_file  = 'ALOS2_P162A_20160720_20161207_Gamma_sim.txt';
% % C band(m) 0.0555041577
% % L band(m) 0.236
% wavelength      = 0.236; 

insar_obs    = load(insar_obs_file);
% geographic postion of the InSAR pixel
lon_obs      = insar_obs(:,1);
lat_obs      = insar_obs(:,2);
% los vector
los_v_obs    = insar_obs(:,3:5);
% los observations
los_obs      = insar_obs(:,6);
unwrap_obs   = -4*pi*los_obs/wavelength;

% local coordinate of InSAR pixel
llh_pixel    = [lat_obs';lon_obs'];
xy_pixel     = llh2localxy(llh_pixel,origin);
x_obs        = xy_pixel(:,1); % East (km)
y_obs        = xy_pixel(:,2); % North (km)

%% 2. Load fault model 
% disp(['...... Load fault model ......']);
% LoadFaultModel;

%% 3. Compute GREENs for RDEs and TDEs
disp(['...... Compute GREENs for RDEs and TDEs ......']);
position_obs = [x_obs,y_obs];
if size(fault_model,2) == 29          % GREENs for RDEs
    G_unproj_obs = compute_RDEs_GREENs(xc,yc,zc,strike,dip,leng,width,position_obs);
elseif size(fault_model,2) == 25      % GREENs for TDEs
    compute_TDEs_GREENs_formula = input('    [Choose the formula to compute TDEs GREENs: 1. OKADA  2. MEADE]');
    if compute_TDEs_GREENs_formula == 1 
        % Using OKADA's formula to compute GREENs function due to TDEs
        G_unproj_obs = compute_TDEs_GREENs_Okada(xc,yc,zc,strike,dip,area,vertices,position_obs);
        % save G_unproj_Okada.txt G_unproj -ascii;
    elseif compute_TDEs_GREENs_formula == 2
        % Using Meade's formular to compute GREENs function due to TDEs
        G_unproj_obs = compute_TDEs_GREENs_Meade(xc,yc,zc,strike,dip,area,vertices,position_obs);
        % save G_unproj_Meade.txt G_unproj -ascii;
    end         
end 

%% 3. Project GREENs of InSAR onto the LOS 
disp('...... Project GREENs of InSAR onto the LOS ......');
% number of the observed InSAR pixel
num_los_obs = numel(los_obs);
G_proj_obs  = zeros(num_los_obs,2*n_patches);
index_tmp   = 1;
for ii=1:num_los_obs
    G_proj_obs(ii,:) = los_v_obs(ii,:) * G_unproj_obs(index_tmp:index_tmp+2,:);
    index_tmp        = index_tmp + 3;
end

%% 4. Simulate the InSAR observations
disp('...... Simulate the InSAR observations ......');
los_sim    = G_proj_obs * slip;
unwrap_sim = -4*pi*los_sim/wavelength;

%% 5. Rewrap the observed and simulated LOS displacements
disp('...... Rewrap the simulated LOS displacements ......');
if (~isreal(unwrap_obs)) helphelp; return; end;
if (~isreal(unwrap_sim)) helphelp; return; end;
wrap_obs = mod(unwrap_obs + pi,2*pi) - pi;
wrap_sim = mod(unwrap_sim + pi,2*pi) - pi;

%% 5. Plotting
disp('...... Plotting ......');
% plot observed, simualted and the residual of LOS displacement
figure;
subplot(1,3,1);
scatter(lon_obs,lat_obs,10,los_obs,'filled');
box on;axis on;axis image;colormap(jet);set(get(colorbar,'ylabel'),'string','LOS(m)');
xlabel('East(km)');ylabel('North(km)');title('Observed InSAR');
caxis([-0.10 0.10]);
subplot(1,3,2);
scatter(lon_obs,lat_obs,10,los_sim,'filled');
box on;axis on;axis image;colormap(jet);set(get(colorbar,'ylabel'),'string','LOS(m)');
xlabel('East(km)');ylabel('North(km)');title('Simulated InSAR');
caxis([-0.10 0.10]);
subplot(1,3,3);
scatter(lon_obs,lat_obs,10,los_obs - los_sim,'filled');
box on;axis on;axis image;colormap(jet);set(get(colorbar,'ylabel'),'string','LOS(m)');
xlabel('East(km)');ylabel('North(km)');title('Residual InSAR');
caxis([-0.10 0.10]);

% plot the observed, simulated and the residual of wrapped interferograms
figure;
subplot(1,3,1);
scatter(lon_obs,lat_obs,10,wrap_obs,'filled');
box on;axis on;axis image;colormap(jet);set(get(colorbar,'ylabel'),'string','LOS(m)');
xlabel('East(km)');ylabel('North(km)');title('Observed InSAR');
subplot(1,3,2);
scatter(lon_obs,lat_obs,10,wrap_sim,'filled');
box on;axis on;axis image;colormap(jet);set(get(colorbar,'ylabel'),'string','LOS(m)');
xlabel('East(km)');ylabel('North(km)');title('Simulated InSAR');
subplot(1,3,3);
scatter(lon_obs,lat_obs,10,wrap_obs - wrap_sim,'filled');
box on;axis on;axis image;colormap(jet);set(get(colorbar,'ylabel'),'string','LOS(m)');
xlabel('East(km)');ylabel('North(km)');title('Residual InSAR');

%% 6. Save the simulated results
disp('...... Save the simulated results ......');
out = [lon_obs,lat_obs,wrap_sim];
save(insar_sim_file,'out','-ascii');





