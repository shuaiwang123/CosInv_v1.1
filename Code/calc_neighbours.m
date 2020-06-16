function Nearest_Points_index = calc_neighbours(data,N_nearest)
% this function is used to compute nearest N_nearest points given locations
% in data, returns the indexes of the nearest N_nearest points to each
% point in the dataset data
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n_points = size(data,1);
%
distance_matrix = zeros(n_points);
% 
% find the distance between every two points, obtain upper triangular
% distance matrix first
for ii=1:n_points-1
    for jj=ii+1:n_points
        distance_matrix(ii,jj) = dist_fnc(data(ii,:),data(jj,:));
    end
end

% fill the lower triangular distance matrix, relying on the fact that
% distance_matrix(ii,ii)=0
distance_matrix = distance_matrix + distance_matrix';

% sort the distances by ascending order using the MATLAB function sort,keep
% the index, and pull out the indexes of the N_nearest nearest points to
% each point of the fault
[distance_temp,distance_index]=sort(distance_matrix,2);
clear('distance_temp');
Nearest_Points_index = distance_index(:,2:1+N_nearest);

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function dist = dist_fnc(x1,x2)
% find the Eucliean norm between two input vectors
dist = norm(x1-x2);