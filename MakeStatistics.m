% This script is used to calculate the SLIP POTENCY, SEISMIC MOMENT(Mo) and
% MOMENT MAGNITUDE(Mw) due to slip at depth as well as the residual between
% observed data and synthetic data. 
%
% The basic flow is as follows
% 1. CALCULATE SEISMIC MOMENT AND SEISMIC MAGNITUDE 
% 2. CALCULATE RESUDUAL BETWEEN OBSERVED AND SYNTHETIC DATA
% 3. MAKEING STATISTICS BETWEEN OBSERVED AND SYNTHETIC DATA
%    S 1. Make statistics about RMS
%    S 2. Make statistics about Residual
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
disp('...... Make Statistics ......');
% [ 1. Calculate seismic moment and seismic magnitude ]
Mo = 0;
Mw = 0;
if size(fault_model,2) == 29
    area    = fault_model(:,6).*fault_model(:,7);
    [Mo,Mw] = cal_Mw(area,Tslip);
elseif size(fault_model,2) == 25
    area    = fault_model(:,6);
    [Mo,Mw] = cal_Mw(area,Tslip);
end

% [ 2. Calculate residual between observed and synthetic data ]
disp('    [Calculate residual between observed and synthetic data ]');
% Residual between observed and synthetic data
res_disp_obs_syn = cell(n_datasets,1);

% GPS residual
res_disp_obs_syn{1} = [gps.disp(:,1) - syngpse,  gps.disp(:,2) - syngpsn,  gps.disp(:,3) - syngpsu];

% InSAR residual
for ii=2:n_datasets
    res_disp_obs_syn{ii} = d(d_index{ii}) - dmodel(d_index{ii});
end

% [ 3. Make statistics between observed and synthetic data ]
% ... S 1. Make statistics about RMS ...
disp('    [Making statistics about the RMS of residual]');

rms_gpse = sqrt(sum((gps.disp(:,1) - syngpse).^2)/gps.nobs);
rms_gpsn = sqrt(sum((gps.disp(:,2) - syngpsn).^2)/gps.nobs);

obs_gps_hori = sqrt(gps.disp(:,1).^2 + gps.disp(:,2).^2);
syn_gps_hori = sqrt(syngpse.^2 + syngpsn.^2);
rms_gps_hori = sqrt(sum((obs_gps_hori - syn_gps_hori).^2)/numel(obs_gps_hori));

rms_gpsu = sqrt(sum((gps.disp(:,3) - syngpsu).^2)/gps.nobs);

rms_insar = cell(n_insar_datasets,1);
for jj=2:n_datasets
    % rms_insar{jj-1} = sqrt(sum((res_disp_obs_syn{jj}).^2)/numel(d_index{jj}));
    rms_insar{jj-1} = rms(res_disp_obs_syn{jj});
end

disp(['      ... The RMS of GPS East components: ',num2str(rms_gpse)]);
disp(['      ... The RMS of GPS North components: ',num2str(rms_gpsn)]);
disp(['      ... The RMS of GPS horizontal components: ',num2str(rms_gps_hori)]);
disp(['      ... The RMS of GPS Up components: ',num2str(rms_gpsu)]);
for ii=1:n_insar_datasets
    disp(['      ... The RMS of InSAR',num2str(ii),': ',num2str(rms_insar{ii})]);
end
% S 1. END

% ... S 2. Make statistics about Residual ...
disp('    [Making statistics about the Residual]');

res_gpse = gps.disp(:,1) - syngpse;
res_gpsn = gps.disp(:,2) - syngpsn;
res_gps_hori = obs_gps_hori - syn_gps_hori;
res_gpsu = gps.disp(:,3) - syngpsu;

res_insar = cell(n_insar_datasets,1);
for jj=2:n_datasets
    res_insar{jj-1} = res_disp_obs_syn{jj};
end

disp(['      ... The max/min/average residual of GPS East components: ',num2str(max(abs(res_gpse))),'  ',num2str(min(abs(res_gpse))),'  ',num2str(mean(res_gpse))]);
disp(['      ... The max/min/average residual of GPS North components: ',num2str(max(abs(res_gpsn))),'  ',num2str(min(abs(res_gpsn))),'  ',num2str(mean(res_gpsn))]);
disp(['      ... The max/min/average residual of GPS horizontal components: ',num2str(max(abs(res_gps_hori))),'  ',num2str(min(abs(res_gps_hori))),'  ',num2str(mean(res_gps_hori))]);
disp(['      ... The max/min/average residual of GPS Up components: ',num2str(max(abs(res_gpsu))),'  ',num2str(min(abs(res_gpsu))),'  ',num2str(mean(res_gpsu))]);
for ii=1:n_insar_datasets
    disp(['      ... The max/min/average residual of InSAR',num2str(ii),': ',num2str(max(abs(res_insar{ii}))),'  ',num2str(min(abs(res_insar{ii}))),'  ',num2str(mean(res_insar{ii}))]);
end






    
    