function G = compute_TDEs_GREENs_Meade(x,y,z,strike,dip,area,vertices,position)
% Computes the Green functions for triangle sources with 6 parameters x,y,z
% strike,dip,area,vertices, and position 
%
% ###
% x : the x(EW) position of center triangles/points,      unit[km]
% y : the y(NS) position of center triangles/points,      unit[km]
% z : the z(Down is positive) position of triangles/points, unit[km]
%
% strike : the strike of triangles/points,    unit[degree]
% dip    : the dip angle of triangles/points, unit[degree]
% length : the length of triangles/points,    unit[km^3]
% width  : the width of triangles/points,     unit[km^3]
%
% position : the position of observation points[in degrees] in local
%               xy coordinates
% ### 
%   if you want to know the format of Green functions computed by the
%   Fortran code, you can refer to ReadMe.txt in GREENFUNC/bin
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
disp('    [Meade formula building TDEs GREENs]');
% define the number of surface displacement directions(3=ENU)
n_surface_disp_directions = 3;
% define the number of fault slip directions(2=SS,DS)
n_slip_directions = 2;

% [ define the observation points on the surface ]
% observation_points = extract_from_cell(position);
observation_points = position;
n_observation_points = size(observation_points,1);
% pad the observation_points matrix with zeros for z=0, as this coordinate
% is not given in position but is required for Okada
observation_points = [observation_points,zeros(n_observation_points,1)];
% the coordinates of observation points
x_obs = observation_points(:,1);
y_obs = observation_points(:,2);
z_obs = observation_points(:,3);

% the number of triangles
n_triangles = size(x,1);

% nu is the Poisson's ratio 
% u/(lamada+u)=1-2*nu
nu     = 0.25;		
% 

G = zeros(n_observation_points*n_surface_disp_directions,...
    n_triangles*n_slip_directions);

% ### Note that only we put a Negative Sign on the right side of equations
%     which marked by(1),(2) and (3) below, can we share the same GREENs
%     function computed by OKADA and MEADE, repectively.
U_0 = 0;
U_1 = -1;     % (1)

% [ compute the GREENs function due to triangle dislocation based on the  ]
% [ formula derived by Meade                                              ]
for ii=1:n_triangles
    % x_tri the x coordinates of triangles vertices
    % y_tri the y coordinates of triangles vertices
    % z_tri the z coordinates of triangles vertices
    x_tri(1) = vertices(ii,1);x_tri(2) = vertices(ii,4);x_tri(3) = vertices(ii,7);
    y_tri(1) = vertices(ii,2);y_tri(2) = vertices(ii,5);y_tri(3) = vertices(ii,8);
    z_tri(1) = vertices(ii,3);z_tri(2) = vertices(ii,6);z_tri(3) = vertices(ii,9);
    %
    Us = CalcTriDisps(x_obs,y_obs,z_obs,x_tri,y_tri,z_tri,nu,U_1,U_0,U_0);
    Ud = CalcTriDisps(x_obs,y_obs,z_obs,x_tri,y_tri,z_tri,nu,U_0,U_0,U_1);
    
    %Build the Green function       
    G(1:n_surface_disp_directions:end,ii)= Us.x;
    G(2:n_surface_disp_directions:end,ii)= Us.y;
    G(3:n_surface_disp_directions:end,ii)= -Us.z;  % (2)
    G(1:n_surface_disp_directions:end,ii+n_triangles)= Ud.x;
    G(2:n_surface_disp_directions:end,ii+n_triangles)= Ud.y;
    G(3:n_surface_disp_directions:end,ii+n_triangles)= -Ud.z;  % (3)
end

