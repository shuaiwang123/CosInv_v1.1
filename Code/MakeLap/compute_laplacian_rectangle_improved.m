function Lap = compute_laplacian_rectangle_improved(x,y,z,vertices)
% This function is used to construct the Laplacian smoothing matrix for
% methods using rectangle dislocations with a impoved scheme
%
% ###
%  The Laplacian operator used in the RDE(rectangle dislocation element)
%  was first put forward by Jonsson et al(2002) in BSSA. The formula of
%  Laplacian operator has been represent by many articles, here you can
%  refer to Jiang et al(2013) in GJI, the [4] formula gives the detail
%  Laplacian operator for RDE. However, generally, there are four adjacent
%  patches for an RDE, but there are only three adjacent patches for the
%  boundary RDEs and two for the RDEs at each of the four corners of the
%  fault.For these corner RDEs, doubling the slip in the numberators of
%  eq(4) in Jiang et al(2013) in GJI means that the slip must be halved to
%  obtain the expected smoothness between two adjacent RDEs, which is
%  obviously unreasonable. For this reason, Jiang et al(2013) modified the
%  classical Laplacian smoothing operator for the boundary RDEs
%
% ### input
%   x        : the x direction[EW] coordinates of the rectangle central
%   y        : the y direction[NS] coordinates of the rectangle central
%   z        : the z direction coordinates of the rectangle central,down is
%              positive
%   vertices : the four vertex coordinates of the rectangle fault model
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
disp('      ... The Improved RDEs Laplacian method has been choosed ...');
%% constructing the vertice index of the rectangles
eps = 1e-3;
n_patches = size(vertices,1);
% t1 : vertex 1 of rectangle patches
% t2 : vertex 2 of rectangle patches
% t3 : vertex 3 of rectangle patches
% t4 : vertex 4 of rectangle patches
index_first_vertex  = 1:3;
index_second_vertex = 4:6;
index_third_vertex  = 7:9;
index_four_vertex   = 10:12;
t1 = vertices(:,index_first_vertex);
t2 = vertices(:,index_second_vertex);
t3 = vertices(:,index_third_vertex);
t4 = vertices(:,index_four_vertex);

t = [t1;t2;t3;t4];
% the unique vertexs on the fault surface
t = unique(t,'rows');
% index for rectangle vertex with coordinates in t
vertex_index = zeros(n_patches,4);
for ii=1:n_patches
    for jj=1:length(t)
        delta_t1 = norm(t1(ii,:)-t(jj,:));
        delta_t2 = norm(t2(ii,:)-t(jj,:));
        delta_t3 = norm(t3(ii,:)-t(jj,:));
        delta_t4 = norm(t4(ii,:)-t(jj,:));
        % determine the index of the first vertex of rectangle with coordinates in t
        if delta_t1 < eps
            vertex_index(ii,1) = jj;
        end
        % determine the index of the second vertex of rectangle with coordinates in t
        if delta_t2 < eps
            vertex_index(ii,2) = jj;
        end
        % determine the index of the third vertex of rectangle with coordinates in t
        if delta_t3 < eps
            vertex_index(ii,3) = jj;
        end
        % determine the index of the four vertex of rectangle with coordinates in t
        if delta_t4 < eps
            vertex_index(ii,4) = jj;
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

%%  constructing the Laplacian
% ### 
%   Note that: For simplicity I assume the distances between adjacent 
%              patches in the row and column directions are equal, and I
%              substitute the distances between adjacent patches in the row
%              and column with the mean value of distances between the
%              interest pathces and its adjacent patches. In general, it is
%              reasonable as the size of the adjacent dislocation is not
%              significantly different.
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
    
    if n_neighbours == 2        % for corner RDEs
        delta_l = mean(distance_matrix(k,common)); % the mean value of distances between the interest pathces and its adjacent patches
        Lap(k,common) =  1./delta_l^2;
        Lap(k,k)      = -2./delta_l^2;
    elseif n_neighbours == 3    % for top and bottom RDEs
         delta_l = mean(distance_matrix(k,common));
         Lap(k,common) =  1./delta_l^2;
         Lap(k,k)      = -3./delta_l^2;
    elseif n_neighbours == 4    % for central RDEs
        delta_l = mean(distance_matrix(k,common));
        Lap(k,common) =  1./delta_l^2;
        Lap(k,k)      = -4./delta_l^2;         
    else
        error('      ... The calculted number of neighbours is not 2 or 3 or 4 in the process of Classical RDEs Laplacian calculted ...');
    end
end
               
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function dist = dist_fnc(x1,x2)
% finds the Eucliean distance between two input parameters
dist = norm(x1-x2);        
        