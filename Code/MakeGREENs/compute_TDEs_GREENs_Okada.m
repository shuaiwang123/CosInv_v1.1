function G = compute_TDEs_GREENs_Okada(x,y,z,strike,dip,area,vertices,position)
% Compute Okada Greens function for point source
%
% x : the x(EW) position of center triangle, unit[km]
% y : the y(NS) position of center triangle, unit[km]
% z : the z(Up is position) position of center triangle, unit[km]
%
% strike : the strike of triangles, unit[degree]
% dip    : the dip angle of triangles, unit[degree]
% area   : the area of triangles, unit[km^3]
% vertices : the xyz postion for the ordering triangle vertexes
%
% position : the position of observation points[in degrees] in local
%               xy coordinates
% By: ShuaiWang
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
disp('    [OKADA formula building TDEs GREENs]');
% change unit(from km to meters)
km2m = 1e3;
x = x*km2m;
y = y*km2m;
z = z*km2m;
area = area*km2m^2;

% define the number of point sources
n_point_sources = numel(x);

% define the number of surface displacement directions(3=ENU)
n_surface_disp_directions = 3;
% define the number of slip displacement directions(2=SS,DS)
n_slip_directions = 2;

% define the observation points on the surface
% observation_points = extract_from_cell(position);
observation_points = position;
% observation_points = position
% number of the observation points
n_observation_points = size(observation_points,1);
% pad the observation_points matrix with zeros for z=0, as this coordinate
% is not given in position but is required for Okada
observation_points = [observation_points,zeros(n_observation_points,1)];
observation_points = observation_points*km2m;


% [define the Greenfunction]
G = zeros(n_surface_disp_directions * n_observation_points,n_point_sources * n_slip_directions);

% [compute the alpha]
% alpha=(lamada+u)/(lamada+2u)=1-(vs/vp)^2
% Vp = 4.5;   % unit:km/s
% Vs = 2.6;   % unit:km/s
% alpha = 1-(Vs/Vp).^2;
GPa = 1e9;
nu  = 0.25;
mu  = 30*GPa;
lambda = 2*mu*nu/(1-2*nu);
alpha  = (lambda+mu)/(lambda+2.0*mu);


% [compute the Greenfunction]
for ii=1:n_point_sources
    % change the refence frame so that the stations are referred to an
    % origin corresponding to middle point of top edge, if we extend fault
    % plane until the surface
    %
    % new origin
    xorigin = x(ii);% - z(ii)/tand(dip(ii)) * cosd(strike(ii));
    yorigin = y(ii);% + z(ii)/tand(dip(ii)) * sind(strike(ii));
    
    % make a rotation around(xorigin, yorigin) so that the strike would be
    % 0, for the Green Function computation purpose of the origin would be
    % the middle point of top edge
    observation_points_x1 = observation_points(:,1) - xorigin;
    observation_points_y1 = observation_points(:,2) - yorigin;
    Rotation1 = [sind(strike(ii)) cosd(strike(ii));-cosd(strike(ii)) sind(strike(ii))];
    temp = Rotation1*[observation_points_x1';observation_points_y1'];
    % the observation points position in fault coordinate system
    observation_points_x2 = temp(1,:)';
    observation_points_y2 = temp(2,:)';
    
    % Green function due to strike slip point sources
    [SSx1,SSy1,SSz] = Okada_DC3D0(alpha,...
                observation_points_x2,observation_points_y2,0,z(ii),dip(ii),area(ii),0,0,0);
    % Green function due to dip slip point sources
    [DSx1,DSy1,DSz] = Okada_DC3D0(alpha,...
                observation_points_x2,observation_points_y2,0,z(ii),dip(ii),0,area(ii),0,0);  
            
     % projection based on the strike
     %Rotation2 = [cosd(strike(ii)) sind(strike(ii));-sind(strike(ii)) cosd(strike(ii))];
     Rotation2 = [sind(strike(ii)) -cosd(strike(ii));cosd(strike(ii)) sind(strike(ii))];
     % strike-slip
     temp = Rotation2*[SSx1';SSy1'];
     SSx  = temp(1,:);
     SSy  = temp(2,:);
     % dip-slip
     temp = Rotation2*[DSx1';DSy1'];
     DSx  = temp(1,:);
     DSy  = temp(2,:);
     
     % [Build the Green funcion]
     % strike-slip
     G(1:n_surface_disp_directions:end,ii) = SSx;
     G(2:n_surface_disp_directions:end,ii) = SSy;
     G(3:n_surface_disp_directions:end,ii) = SSz;
     % dip-slip
     G(1:n_surface_disp_directions:end,ii+n_point_sources) = DSx;
     G(2:n_surface_disp_directions:end,ii+n_point_sources) = DSy;
     G(3:n_surface_disp_directions:end,ii+n_point_sources) = DSz;
end
     
     
     
     
     
     
     
     
    





















