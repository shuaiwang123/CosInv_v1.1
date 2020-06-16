% function Lap = compute_laplacian_v1(x,y,z,N_nearest)
function Lap = compute_laplacian_irregular(x,y,z,N_nearest)
% this function is used to generate Laplacian with respect to a set of
% points, this function takes in the coordinates of a set of points that
% are randomly scattered on an unknow surface and generates a discrete
% approximation of the Laplacian using the nearest N points.
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
disp('      ... The Irregular Laplacian method for RDEs and TDEs has been choosen ...');
% define the number of points on the surface
n_points = numel(x);
% define the data matrix,where each column is a different coordinate
data = [x,y,z];
% pre-allocate the distance matrix and the Laplacian matrix
Lap = zeros(n_points);
Lap_unscaled = zeros(n_points);

% calculate the indexes of the N_nearest nearest points given locations
% defining a smoothing surface
Nearest_Points_index = calc_neighbours(data,N_nearest);

%% �ޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡ�
% *          Calculate the Laplacian individually for each point         *
% *                                                                      *
% * see:Geertjan Huiskamp. Different formulas for the surface laplacian  *
% * on a triangulated surface. Journal of computational physics, 95(2):  *
% * 477-496,1991. OR the PCAIM manual for details                        *
% �ޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡ�
for ii=1:n_points
    % center the data at data(ii)
    translated_data = data(Nearest_Points_index(ii,:),:) - repmat(data(ii,:),N_nearest,1);
    
    % find the best fitting plane going through data(ii), namely that
    % spanned by the columns of v, and the distance of the points'
    % projection on this plane distances along these axes are u*s
    [u,s,v]=svds(translated_data,2);
    
    % ### actually this will not happen
    % check that size(u,2). if not, recalculate svds(tranlated_data)
    % considering more components
    i_svds = 2;
    while size(u,2) == 1
        [u,s,v]=svds(translated_data,i_svds);
        i_svds = i_svds+1;
    end 
    % ###
    
    % convert to polar coordinates to sort by angle
    %  theta : is the angle with the X axis positive direction,couterclockwise
    %  theta and r are the elements defining the polar coordiantes
    [theta,r] = cart2pol(u(:,1)*s(1,1),u(:,2)*s(2,2));
    %
    % sort by angle
    [theta,theta_index] = sort(theta);
    
    % recorder the points for Nearest_Points_index, r and u
    Nearest_Points_index(ii,:) = Nearest_Points_index(ii,theta_index);
    r = r(theta_index);
    u = u(theta_index,:);
    
    % calculate the differences in theta for use in the Laplacian approximation
    delta_theta = diff([theta;theta(1)]);
    % assign thetas in the positive rotation direction, as the angle from
    % ri to ri+1
    theta_plus = delta_theta;
    % assign thetas in the negative rotation direction, as the angle from
    % ri-1 to ri
    theta_minus = [delta_theta(end);delta_theta(1:end-1)];
    % assign the average distance to the neighbouring points
    r_bar = mean(r);
    
    % calculate theta_tot as in Huiskamp,1991,for the detail formulas that
    % calculate the Theta_tot you can see PCAIM manual P212 D.11 and D.12
    Theta_tot = calc_Theta_tot(theta_minus,theta_plus);
    
    % laplacian approximation, for detail formulas that calculate the
    % laplacian approximation you can see PCAIM manual P212 D.11 and D.12
    for jj=1:N_nearest
        Lap(ii,Nearest_Points_index(ii,jj)) = 4/r_bar*1/Theta_tot*1/r(jj)*...
            calc_Theta(theta_plus(jj),theta_minus(jj));
        
        Lap_unscaled(ii,Nearest_Points_index(ii,jj)) = 4/r_bar*1/Theta_tot*...
            calc_Theta(theta_plus(jj),theta_minus(jj));
    end
    % central point is equal to the negative sum of the neighbouring weights
    Lap(ii,ii) = -sum(Lap(ii,Nearest_Points_index(ii,:)));
    Lap_unscaled(ii,ii) = -sum(Lap_unscaled(ii,Nearest_Points_index(ii,:)));
end

%% �ޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡ�
% * everything below here is new as of MAY 20,2010. aims to fix the      *
% * laplacian corners problem on rectangular faults                      *
% �ޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡޡ�
lap_diag = diag(Lap_unscaled);
lap_med  = median(lap_diag);
lap_std  = std(lap_diag);
% find the bad laplacian index
lap_bad_index = find(abs(lap_diag-lap_med) > 3*lap_std);
for ii=1:numel(lap_bad_index)
    curr_index = lap_bad_index(ii);
    for jj=1:N_nearest
    % the original code
    %    Lap(ii,Nearest_Points_index(ii,jj)) = 4/r_bar*1/Theta_tot*1/r(jj)*...
    %        calc_Theta(theta_plus(jj),theta_minus(jj))*(lap_med/lap_diag(curr_index));
    % the code modified by ShuaiW
       Lap(curr_index,Nearest_Points_index(curr_index,jj)) = 4/r_bar*1/Theta_tot*1/r(jj)*...
           calc_Theta(theta_plus(jj),theta_minus(jj))*(lap_med/lap_diag(curr_index));
    end
    % the original code
    %Lap(ii,ii) = -sum(Lap(ii,Nearest_Points_index(ii,:)));
    % the code modified by ShuaiW
    Lap(curr_index,curr_index) = -sum(Lap(curr_index,Nearest_Points_index(curr_index,:)));
end
    
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function Theta_tot = calc_Theta_tot(theta_plus,theta_minus)
% cal_theta_tot defined as in [Huiskamp,1991]
Theta_tot = 0;
for i=1:numel(theta_plus)
    Theta_tot = Theta_tot + calc_Theta(theta_plus(i),theta_minus(i));
end

function Theta = calc_Theta(theta_plus,theta_minus)
% cal_Theta defined as in [Huiskamp,1991]
Theta = calc_half_Theta(theta_plus) + calc_half_Theta(theta_minus);

function partial_Theta = calc_half_Theta(theta_diff)
% calc_half_Theta defined in [Huiskamp,1991] for theta_plus/theta_minus
tol = 10^(-10);
if (abs(theta_diff) < tol)
    partial_Theta = 0;
else
    partial_Theta = (1-cos(theta_diff))/sin(theta_diff);
end


    
    
    
    
    
    
    
    
    
    
    
    
    
    