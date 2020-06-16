% This script is used to compute GREENs function due to dislocation at
% depth.
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
disp('...... Compute GREENs Due to Dislocation at Depth ......');

%% 1. Compute GREENs for RDEs and TDEs
position = [x,y];
if size(fault_model,2) == 29          % GREENs for RDEs
    G_unproj = compute_RDEs_GREENs(xc,yc,zc,strike,dip,leng,width,position);
elseif size(fault_model,2) == 25      % GREENs for TDEs
    compute_TDEs_GREENs_formula = input('    [Choose the formula to compute TDEs GREENs: 1. OKADA  2. MEADE]');
    if compute_TDEs_GREENs_formula == 1 
        % Using OKADA's formula to compute GREENs function due to TDEs
        G_unproj = compute_TDEs_GREENs_Okada(xc,yc,zc,strike,dip,area,vertices,position);
        % save G_unproj_Okada.txt G_unproj -ascii;
    elseif compute_TDEs_GREENs_formula == 2
        % Using Meade's formular to compute GREENs function due to TDEs
        G_unproj = compute_TDEs_GREENs_Meade(xc,yc,zc,strike,dip,area,vertices,position);
        % save G_unproj_Meade.txt G_unproj -ascii;
    end         
end

%% 2. Project GREENs of InSAR onto the LOS
disp('    [Project GREENs of InSAR onto the LOS]');
n_patches = numel(xc);
G_proj = zeros(numel(d),2*n_patches);

% GREENs of GPS
G_proj(d_index{1},:) = G_unproj(G_unproj_index{1},:);
% GREENs of InSAR
for ii=2:n_datasets
    index = 1;
    for jj=1:numel(d_index{ii})
        G_proj(d_index{ii}(jj),:) = los_v{ii}(jj,:) * G_unproj(G_unproj_index{ii}(index:index+2),:);
        index = index + 3;
    end
end


