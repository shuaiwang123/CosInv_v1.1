% This script is used to compute Laplacian operator for smoothing slip
% gradient. 
% Laplation options
%   'n_neighbours',n_neighbors
%       this option tells the Laplacian how many neightboring points(an
%       integer) to use in its approximation of the Laplician. we recommend
%       somewhere between 4 and 10. try various numbers and see how the
%       results change.
%    'no_slip_points',no_slip_points
%        this option specifies a number of fault elements(integers in a
%        vector) that have heavy penalties applied to any slip on them,
%        this is useful for forcing slip to go to zero at the fault
%        boundaries, it is note that this options overrides
%        'free_surface_depth'
%     'free_surface_depth',free_surface_depth. 
%        1. for any patch whose center is above free_surface_depth, there
%           are no penalties on slip and the standard Laplacian is used.
%        2. for any patch who center is below free_surface_depth, slip is
%           heavily penalized.
%        Note: this option is overridden by 'no_slip_points'
%     'projected'
%        this string should be included if the user wants the algorithm to
%        guess the edge patches(on which slip should be reduced to zero) on
%        the fault using a slight modification of the MATLAB convex hull
%        algorithm
%     'scaling_edge_factor',scaling_edge_factor
%        scaling_edge_factor is the multiplicative factor(doulbe) of a
%        standard distance between patches that is used to guess whether a
%        given patch is on the edge or not.
% ### Note: there are many real earthquakes who can rupture to the crust
%           surface, therefore it is not reasonable to impose the slip on
%           the edge of fault to zero. Furthermore, the high precise of
%           observe datas can given a high resolution of fault rupture and
%           its distribution. So, we can now inverse the fault slip without
%           adding constrains such as imposed zero slip patches .
%           Therefore, we can set the final IND_NO_SLIP_PATCHES to zeros in
%           the following.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
disp('...... Compute Smoothing Operator ......');
smooth_method = input('    [Choose the Smoothing method: 1. Laplacian smoothing  2. Minimum moment]');

%% Find patches lie in the edge, and return their indexes
% === synthetic earthqauke ===
% Give the fault average strike and dip angle which will be used in the
% TEST_FIND_EDGE 
% strike_average = 295;
% dip_average    = 10;
fault_file = load(fault_model_file);

% EDGE patches with no slip
% test_find_edges;
% ind_no_slip_patches = ind_patches_edge';
ind_no_slip_patches = []';

% DEPTH with no slip
% for curved rec
% ind_no_slip_patches_depth = []';
ind_no_slip_patches_depth=find(fault_file(:,3)<=11.09 | fault_file(:,3)>=27.65);

% CENTER patches within the fault surface with no slip
% % for curved rec
ind_center = unique([]');

% FINAL patches with no slip
ind_no_slip_patches = unique([ind_no_slip_patches;ind_no_slip_patches_depth]);
ind_no_slip_patches = unique([ind_no_slip_patches;ind_center]);

%% Laplacian options
laplacian_options = {'n_neighbours',4,'free_surface_depth',100,'projected',...
    'scaling_edge_factor',1,'no_slip_points',ind_no_slip_patches};

%% Get Laplacian operator
if smooth_method == 1
    disp('      ... The Laplacian Operator has been choosed ...');
    [Lap,iedge,laplacian_method_flag] = compute_laplacian_deriver(fault_model,laplacian_options);
elseif smooth_method == 2
    disp('    [Using Minimum moment to smooth slip]');
    Lap = eye(2*n_patches);
end


