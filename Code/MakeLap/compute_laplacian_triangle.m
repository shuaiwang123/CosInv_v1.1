function Lap = compute_laplacian_triangle(vertices)
% This function is used to construct the Laplacian smoothing matrix for
% methods using triangular dislocations.
% 
% ###
% The principle is described in Barnhart et al(2010) in the JGR, in order
% to construct the Laplacian smoothing matrix for methods using triangular
% dislocations, we first identify if the dislocation of interest is in
% contact with either Two(a corner dislocation) or Three(an internal or
% side dislocation) other dislocations. We assign a weight of 1 to the
% dislocation of interest and equivalent values to each adjacent
% dislocation such that the sum of the weight of all dislocations is 0. We
% do no weight each dislocation by its area because we assume the area of
% adjacent dislocations is not significantly different since dislocation
% size should vary smoothly.
% 
%
% ### input
%   vertices : the three vertex coordinates of the triangle fault model
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
disp('      ... The Classical TDEs Laplacian method has been choosed ...');
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
Lap = eye(length(vertex_index));
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
    
    if n_neighbours == 1
        Lap(k,common) = -1;
    elseif n_neighbours == 2
        Lap(k,common) = -1/2;
    elseif n_neighbours ==3
        Lap(k,common) = -1/3;
    else
        error('      ... The calculted number of neighbours is not 1 or 2 or 3 in the process of Classical TDEs Laplacian calculted ...');
    end
end

    
        


    
    
    
    
    
    

