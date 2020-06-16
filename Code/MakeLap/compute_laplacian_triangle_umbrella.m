function Lap = compute_laplacian_triangle_umbrella(x,y,z,vertices)
% This function is used to construct the Laplacian smoothing matrix for
% methods using triangular dislocations.
% 
% ###
% The principle is described in Barnhart et al(2010) in the JGR, in order
% to construct the Laplacian smoothing matrix for methods using triangular
% dislocations, we first identify if the dislocation of interest is in
% contact with either Two(a corner dislocation) or Three(an internal or
% side dislocation) other dislocations. In the classical method to build 
% Laplacian operator for triangle dislocation, a weight of 1 is assigned to 
% the dislocation of interest and equivalent values to each adjacent
% dislocation such that the sum of the weight of all dislocations is 0. We
% do no weight each dislocation by its area because we assume the area of
% adjacent dislocations is not significantly different since dislocation
% size should vary smoothly. However, it is not suite for the model build
% by CFMM as the size of triangle elements varying at depth, instead of a
% finite-difference formulation or classical triangle formulation we use
% the scale-dependent umbrella operator(Desbrun et al(1999), Maerten et
% al(2005) in BSSA).
% 
%
% ### input
%   x        : the x direction[EW] coordinates of the rectangle central
%   y        : the y direction[NS] coordinates of the rectangle central
%   z        : the z direction coordinates of the rectangle central,down is
%              positive
%   vertices : the four vertex coordinates of the rectangle fault model
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
disp('      ... The Scale-base Umbrella TDEs Laplacian method has been choosed ...');
%% constructing the vertice index of the triangles
eps = 1e-3;
n_patches = size(vertices,1);
% t1 : vertex 1 of triangular patches
% t2 : vertex 2 of triangular patches
% t3 : vertex 3 of triangular patches
index_first_vertex  = 1:3; 
index_second_vertex = 4:6;
index_third_vertex  = 7:9;
t1 = vertices(:,index_first_vertex);
t2 = vertices(:,index_second_vertex);
t3 = vertices(:,index_third_vertex);

t = [t1;t2;t3];
% the unique vertexs on the fault surface
t = unique(t,'rows'); 
% index for triangles vertex with coordinates in t
vertex_index = zeros(n_patches,3);
for ii=1:n_patches
    for jj=1:length(t)
        delta_t1 = norm(t1(ii,:)-t(jj,:));
        delta_t2 = norm(t2(ii,:)-t(jj,:));
        delta_t3 = norm(t3(ii,:)-t(jj,:));
        % determine the index of the first vertex of triangle  with coordinates in t 
        if delta_t1 < eps
            vertex_index(ii,1) = jj;
        end
        % determine the index of the second vertex of triangle with coordinates in t
        if delta_t2 < eps
            vertex_index(ii,2) = jj;
        end
        % determine the index of the third vertex of triangle with coordinates in t 
        if delta_t3 < eps
            vertex_index(ii,3) = jj;
        end
    end
end

%% constructing the Laplacian
data = [x,y,z];
distance_matrix = zeros(n_patches);
for k = 1:n_patches-1
    for j = k+1:n_patches
        distance_matrix(k,j) = dist_fnc(data(k,:),data(j,:));
    end
end
distance_matrix = distance_matrix + distance_matrix';

%% constructing the Laplacian 
Lap = zeros(length(vertex_index));
for k=1:length(vertex_index)
    clear common
    common = [];
    a = vertex_index(k,:);
    for j=1:length(vertex_index)
        b = vertex_index(j,:);
        temp = intersect(a,b);
        if length(temp) == 2
            common = [common j];
        end
    end
    n_neighbours = length(common);
    
    L = sum(distance_matrix(k,common)); % the sum of the element center distances, Maerten et al(2005) (4) in BSSA
    
    if n_neighbours == 1 || n_neighbours == 2 || n_neighbours == 3
        Lap(k,common) =  2./(L*distance_matrix(k,common));
        %
        temp = sum(1./distance_matrix(k,common)); %
        Lap(k,k)      = -2/L*temp;
    else
        error('      ... The calculted number of neighbours is not 1 or 2 or 3 in the process of Scale-base Umbrella Laplacian calculted ...');
    end
end

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function dist = dist_fnc(x1,x2)
% finds the Eucliean distance between two input parameters
dist = norm(x1-x2);        

    
        


    
    
    
    
    
    

