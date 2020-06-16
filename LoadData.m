% This script is used to load GPS and InSAR data. Multi-interferograms 
% cover the rupture area can bee loaded by this script. Note that  both 
% observed displacements and variance of GPS and InSAR are in [m] unit. We
% assume there only exist ONE GPS dataset.
%
% GPS.mat file contains a struct GPS and this struct has 7 variables
%  NOBS number of GPS stations
%  LON  longitude of GPS stations
%  LAT  latitude of GPS stations
%  X    East position of GPS stations in local coordinates
%  Y    North position of GPS stations in local coordinates
%  DISP 3D displacements of GPS stations, is a [N*3] matrix, where N is the
%       number of GPS stations. Its format is 
%    east    north    up
%  VAR  variance of GPS observations, is a [N*N] matrix
%
% INSAR.mat file contains a struct INSAR and this struct has 8 variables
%  NOBS number of GPS stations
%  LON  longitude of GPS stations
%  LAT  latitude of GPS stations
%  X    East position of GPS stations in local coordinates
%  Y    North position of GPS stations in local coordinates
%  LOS  displacement along Line-of-Sight, is a [M*1] matrix, where M is the
%       number of InSAR pixels
%  VECT projection vector of LOS displacement, is a [M*3] matrix.
%  VAR  variance of InSAR observations, is a [M*M] matrix.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% global n_datasets d_index position_index lon lat x y los_v
% global d covd rampg nramp
%%
disp('...... Load Observed Data ......');
% Number of GPS datasets
n_gps_datasets = 1;
% Number of interferograms
n_insar_datasets = 2;
% Number of datasets including GPS and InSAR
n_datasets = n_gps_datasets + n_insar_datasets;
% Rescaler to determine the relative weight between each datasets
rescale = {1,1,1,1,1}; % check board

d    = [];
covd = [];
lon  = [];
lat  = [];
x    = [];            % E-W[km]
y    = [];
los_v = cell(n_datasets,1);
% D_INDEX is an index pointing towards the indices of d relative to the
% observation points of GPS and InSAR 
d_index = cell(n_datasets,1);
% POSITION_INDEX is an index pointing towards the indices of lon,lat,x and
% y relative to the postion of GPS stations and InSAR pixels
position_index = cell(n_datasets,1);
% G_UNPROJ_INDEX is an index pointing towards the indices of unprojected GREENs 
% relate to the GPS stations and InSAR pixels
G_unproj_index = cell(n_datasets,1);

%% 1. Load GPS
if n_gps_datasets ~= 0
    disp('    [Load GPS Data]');
    load('synthetic_GPS.mat');
    
    % Prelocate the GPS observations to row vector
    d          = [d; reshape(gps.disp',gps.nobs*3,1)];
    covd       = blkdiag(covd,gps.var/rescale{1});
    los_v{1}   = [];
    d_index{1} = 1:gps.nobs*3;
    G_unproj_index{1} = 1:gps.nobs*3;
    
    lon        = [lon;gps.lon];
    lat        = [lat;gps.lat];
    x          = [x;gps.x];
    y          = [y;gps.y];
    position_index{1} = 1:gps.nobs;
end

%% 2. Load InSAR
if n_insar_datasets ~= 0 
    for ii=1:n_insar_datasets
        disp(['    [Load InSAR',num2str(ii),' Data]']);
        load(['synthetic_InSAR',num2str(ii)]);
        
        d             = [d; insar.los];
        %insar.var     = eye(insar.nobs);
        covd          = blkdiag(covd,insar.var/rescale{ii+1});
        los_v{ii+1}   = insar.vect;
        d_index{ii+1} = d_index{ii}(end)+1:d_index{ii}(end)+insar.nobs;
        G_unproj_index{ii+1} = G_unproj_index{ii}(end)+1:G_unproj_index{ii}(end)+insar.nobs*3;
        
        lon = [lon;insar.lon];
        lat = [lat;insar.lat];
        x   = [x;insar.x];
        y   = [y;insar.y];
        position_index{ii+1} = position_index{ii}(end)+1:position_index{ii}(end)+insar.nobs;
    end
end

%% 3. Create ramp matrix for InSAR
disp('...... Create Ramp Matrix for InSAR ......');
% ramp_shape = input('    [Choose the shape of ramps to be inverted for : 1. No ramp  2. Linear ramp  3. Bilinear ramp  4. Quadratic ramp]'); 
ramp_shape = 1;
% Ramp parameters
rampg =[];
% If there exist InSAR observed data
if n_datasets >=2
    for ii=2:n_datasets
        id    = position_index{ii};
        % Xtmp  = x(id) - min(x(id));
        % Ytmp  = y(id) - min(y(id));
        Xtmp  = x(id) - mean(x(id));
        Ytmp  = y(id) - mean(y(id));
        XXtmp = Xtmp.^2;
        YYtmp = Ytmp.^2;
        XYtmp = Xtmp.*Ytmp;
        XXtmp = XXtmp/max(XXtmp);
        YYtmp = YYtmp/max(YYtmp);
        XYtmp = XYtmp/max(XYtmp);
       %% Comment and uncomment the different options below based on the shape of ramps to be inverted for
        if ramp_shape == 1,disp('    [No ramp]');
            rampg = []; % No ramp
        elseif ramp_shape == 2,disp('    [Linear ramp]');
            rampg = blkdiag(rampg,[Xtmp Ytmp ones(numel(id),1)]); % Linear ramp
        elseif ramp_shape == 3,disp('    [Bilinear ramp]');
            rampg = blkdiag(rampg,[Xtmp Ytmp ones(numel(id),1) XYtmp XXtmp YYtmp]); % Bilinear ramp
        else disp('    [Quadratic ramp]');
            rampg = blkdiag(rampg,[Xtmp Ytmp ones(numel(id),1) Xtmp.*Ytmp Xtmp.^2 .*Ytmp.^2]); % Quadratic ramp
        end
    end
    % For GPS we fill RAMPG with 0 to keep the demension of RAMPG consistent
    % with GREENs
    nramp = size(rampg,2);
    rampg = [zeros(gps.nobs*3,nramp);rampg];
else
% If there is only GPS observed data
    nramp = size(rampg,2);
end


    
    
    