function [Lap,iedge,laplacian_method_flag] = compute_laplacian_deriver(fault_model,options)
% This function takes in a fault model, computes an approximation of the
% Laplacian for the fault model, and attempt to find the edges of the fault
% and assigns these to [iedge].
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
disp('    [Using Laplacian operator to smooth slip]');
n_options = numel(options);
kk = 1;

%% default options
projection_flag        = 0;
edge_points_given_flag = 0;
zmax                   = NaN;
scaling_edge_factor    = 1;
strike                 = fault_model(:,4);
dip                    = fault_model(:,5);

%% load fault model
xc = fault_model(:,1);
yc = fault_model(:,2);
zc = fault_model(:,3);
% the number of sub-fault patches[can be point,rectangle and triangle
% elements]
n_patches = numel(xc);

%% load options
while kk<=n_options
    switch options{kk}
        case 'n_neighbours' % the number of neighbours to use for the Laplacian approximation
            n_neighbours = options{kk+1};
            if ~isnumeric(n_neighbours)
                error('    [number of neighbours not given]');
            end
            kk = kk+1;
        case 'no_slip_points' % slip will be penalized on those points given in an index or file
            no_slip_patches = options{kk+1};
            if ischar(no_slip_patches)
                iedge = load(no_slip_patches);
                if size(iedge,2)>size(iedge,1), iedge=iedge';end
            else
                iedge = no_slip_patches(:);
            end
            edge_points_given_flag = 1;
            kk = kk+1;
        case 'free_surface_depth' % depth shallower than which slip is allowed because we have a free surface
            zmax=options{kk+1};
            if ~isnumeric(zmax)
                error('    [limit depth of free surface not given]');
            end
            kk = kk+1;
        case 'projected' % do you want to project onto the best-fitting plane to guess the edges=
            projection_flag = 1;
        case 'scaling_edge_factor' % multiplicative factor describing how close to the edge of convex huxx a dislocation must be to count as an edge point
            scaling_edge_factor = options{kk+1};
            if ~isnumeric(scaling_edge_factor)
                error('    [scaling edge factor not given]');
            end
            kk = kk+1;
%         case 'strike_angle' % strike angle for each patch
%             strike = options{kk+1};
%             if ~isnumeric(strike)
%                 error('    [strike angle not given]');
%             end
%             kk = kk+1;
%         case 'dip_angle' % dip angle for each patch
%             dip = options{kk+1};
%             if ~isnumeric(dip)
%                 error('    [dip angle not given]');
%             end
%             kk = kk+1;           
    end
    kk = kk+1; % increment to next value
end

%% calculation of trigonometrics parameters
median_strike = median(strike);
median_dip    = median(dip);
sin_str       = sind(median_strike);
cos_str       = cosd(median_strike);
sin_dip       = sind(median_dip);
cos_dip       = cosd(median_dip);

%% if the edge patches on which the slip will be heavily penalties, 
%% we will use the Convex hull function in MATLAB to find the edge
if edge_points_given_flag == 0
    % [ compute distances between dislocations ]
    distance_matrix = zeros(n_patches);
    data = [xc,yc];
    for ii = 1:n_patches-1
        for jj = ii+1:n_patches
            distance_matrix(ii,jj) = dist_fnc(data(ii,:),data(jj,:));
        end
    end
    distance_matrix = distance_matrix + distance_matrix';
    sorted_distance = sort(distance_matrix,2);
    dl0 = mean(sorted_distance(:,2));% estimate of the mean distance between points
    dl  = scaling_edge_factor*dl0; % assume that all points within +/- fact_edge_out*dl/2 are edge points
    %######################################################################
    if projection_flag ==1
        % [ project fault on mean local axes ]
        u1 = zeros(3,1);
        u2 = zeros(3,1);
        u3 = zeros(3,1);

        u1(1) = sin_str;
        u1(2) = cos_str;
        u1(3) = 0;
        
        u2(1) =  cos_str*cos_dip;
        u2(2) = -sin_str*cos_dip;
        u2(3) =  sin_dip;
        
        u3(1) =  cos_str*sin_dip;
        u3(2) = -sin_str*sin_dip;
        u3(3) = -cos_dip;
        
        x1 = xc.*u1(1) + yc.*u1(2) + zc.*u1(2);
        y1 = xc.*u2(1) + yc.*u2(2) + zc.*u2(2);
        
        icontour = convhull(x1,y1); % use the convex hull algorithm to find the edge patches
        x2 = x1(icontour); % the edge patches derived by convex hull algorithm
        y2 = y1(icontour); % the edge patches derived by convex hull algorithm
        n1 = length(x1); % the number of patches
        n2 = length(x2); % the number of edge patches derived by Convec hull algorithm
        
        % [ compute distance between the edges given by Convexhull and    ]
        % [ pick the points closer than [dl] pick                         ]
        M1 = [x2(1),y2(1)];
        junk_index = [];
        % To find the edge patches based on the convex hull points derived
        % by convex hull algorithm
        for ii=2:n2
            M2 = [x2(ii),y2(ii)];
            dM = M2-M1;
            dM = dM/norm(dM,'fro');
            M1 = M2;
            distance_vector = zeros(n1,1);
            for jj=1:n1
                Mcurr = [x1(jj),y1(jj)];
                McurrM1 = Mcurr-M1;
                distance_vector(jj) = sqrt(McurrM1*McurrM1' - dot(McurrM1,dM)^2);
            end
            igood = find(distance_vector<dl/2); % find the patches which lies in the edge
            junk_index = [junk_index;igood];
            %plot(x1(igood),y1(igood),'ms');
        end
        iedge = unique(junk_index);
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    elseif projection_flag == 0
        % [ find the edge assuming a rectangle fault with dip angle>0 ]
        [x1,y1] = rotate2d_center(xc,yc,mean(xc),mean(yc),median_strike);
        x1_min = min(x1);
        x1_max = max(x1);
        y1_min = min(y1);
        y1_max = max(y1);
        iedge = find(abs(x1-x1_min)<dl/2 | abs(x1-x1_max)<dl/2 | abs(y1-y1_min)<dl/2 | abs(y1-y1_max)<dl/2);
    end
    %######################################################################
end

%% ========================================================================
%   *           exclude from index files points shallow than zmax         *
%  ==============================================================================
if ~isnan(zmax)
    ifree = find(zc<=zmax);
    if edge_points_given_flag == 1
        no_free = intersect(ifree,iedge);
        ifree   = setdiff(ifree,no_free);
    end
    iedge = setdiff(iedge,ifree);
end

%% ========================================================================
%   *  first compute the Laplacian and then excluding points in index edge
%   ============================================================================
disp('      ... 1. Classical RDEs  2. Improved RDEs  3. Classical TDEs  4.Scale-dependent umbrella TDEs  5.Irregular both for RDEs and TDEs ...');
laplacian_method_flag = input('      ... Please choose the method to compute the Laplacian Operator ...');

vertices_index_in_rec_fault_model = 8:19;
vertices_index_in_trg_fault_model = 7:18;
if laplacian_method_flag == 1
    vertices = fault_model(:,vertices_index_in_rec_fault_model);
    if exist('Lap_matrix.mat','file')
        load('Lap_matrix.mat');
        Lap = Lap;
    elseif ~exist('Lap_matrix.mat','file')
        Lap = compute_laplacian_rectangle(xc,yc,zc,vertices);
        save('Lap_matrix.mat','Lap');
    end
elseif laplacian_method_flag == 2
    vertices = fault_model(:,vertices_index_in_rec_fault_model);
    if exist('Lap_matrix.mat','file')
        load('Lap_matrix.mat');
        Lap = Lap;
    elseif ~exist('Lap_matrix.mat','file')
        Lap = compute_laplacian_rectangle_improved(xc,yc,zc,vertices);
        save('Lap_matrix.mat','Lap');
    end
elseif laplacian_method_flag == 3
    vertices = fault_model(:,vertices_index_in_trg_fault_model);
    if exist('Lap_matrix.mat','file')
        load('Lap_matrix.mat');
        Lap = Lap;
    elseif ~exist('Lap_matrix.mat','file')
        Lap = compute_laplacian_triangle(vertices);
        save('Lap_matrix.mat','Lap');
    end    
elseif laplacian_method_flag == 4
    vertices = fault_model(:,vertices_index_in_trg_fault_model);
    if exist('Lap_matrix.mat','file')
        load('Lap_matrix.mat');
        Lap = Lap;
    elseif ~exist('Lap_matrix.mat','file')
        Lap = compute_laplacian_triangle_umbrella(xc,yc,zc,vertices);
        save('Lap_matrix.mat','Lap');
    end    
elseif laplacian_method_flag == 5
    vertices = fault_model(:,vertices_index_in_trg_fault_model);
    if exist('Lap_matrix.mat','file')
        load('Lap_matrix.mat');
        Lap = Lap;
    elseif ~exist('Lap_matrix.mat','file')
        Lap = compute_laplacian_irregular(xc,yc,zc,n_neighbours);
        save('Lap_matrix.mat','Lap');
    end    
end

lag = median(diag(Lap));

% including strike_slip Laplacian smoothing and dip-slip Laplacian smoothing
Lap = [Lap,zeros(size(Lap));zeros(size(Lap)),Lap]; 

Lap_noslip_edge = Lap;
n_patches_edge  = length(iedge);

for ii=1:n_patches_edge
    index1 = iedge(ii);          % for strike-slip
    index2 = index1 + n_patches; % for dip-slip
    % ###
    %   Note that some earthquakes can rupture the surface, so there is no
    %   needs to do the following
    Lap_noslip_edge(index1,:) = 0;
    Lap_noslip_edge(index2,:) = 0;
    % ###
    %   proper way to penalize slip, But may not work if the number of
    %   neighbours is not properly chosen. [See trike below].
%     Lap_noslip_edge(index1,index1) = Lap(index1,index1);
%     Lap_noslip_edge(index2,index2) = Lap(index2,index2);   
%
    % ### trike to fix the problem of a bad choice of neighbours
    Lap_noslip_edge(index1,index1) = lag;
    Lap_noslip_edge(index2,index2) = lag;
end
Lap = Lap_noslip_edge;

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function dist = dist_fnc(x1,x2)
% finds the Eucliean distance between two input parameters
dist = norm(x1-x2);
    
        
        
            

        












