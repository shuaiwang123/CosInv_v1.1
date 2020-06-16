function G = compute_RDEs_GREENs(x,y,z,strike,dip,length,width,position)
%COMPUTE_MATLAB_RECTANGULAR_SOURCE   Compute Okada Greens function for
%                                    rectangular patches
%
%   G = COMPUTE_MATLAB_RECTANGULAR_SOURCE(X,Y,Z,STRIKE,DIP,LENGTH,WIDTH,
%   POSITION) computes the Greens functions for 
%   patches with 7 parameters X, Y, Z, STRIKE, DIP, LENGTH, WIDTH and 
%   positions on the surface given in POSITION. 
%
%   Example:
%   PCAIM_driver
%
%   See also get_fault_model,compute_point_source, PCAIM_driver.

%   By Marion Thomas
%   Copyright 2009-2010 Tectonics Observatory
%   $Revision: 1.0.0.0 $  $Date: 2010/06/29  $

disp('    [OKADA formula building RDEs GREENs]');
%Change UNIT (from km to meters)
x=x*1e3;y=y*1e3;z=z*1e3;length =length*1e3;width = width*1e3;

% Define the number of patches
n_patches = numel(x);

% Define the number of surface displacement directions (3 = ENU)
n_surface_disp_dirs = 3;
% Define the number of slip displacement directions (2 = SS,DS)
n_slip_directions = 2;

% Define the observation points on the surface
% observation_points = extract_from_cell(position);
observation_points = position;
% observation_points =position;
n_observation_points = size(observation_points,1);
% Pad the observation_points matrix with zeros for z = 0, as this
% coordinate is not given in position but is required for Okada
observation_points = [observation_points, zeros(n_observation_points,1)];


%Define variables
G=zeros(n_surface_disp_dirs* n_observation_points,...
    n_patches * n_slip_directions);
% u/(lamada+u)=1-2*nu
nu     = 0.25;		%Poisson's ratio 

%Compute Green Function

for k=1:n_patches
    
    %Change the refence frame so that the stations are referred to an origin
    %corresponding to middle point of top edge, if we extend fault plane 
    %until the surface
    
    %new origin
    xorigin = x(k) -z(k)/tand(dip(k))* cosd(strike(k));
    yorigin = y(k) +z(k)/tand(dip(k))* sind(strike(k));
    
    %Make a roation around (xorigin,yorigin) so that the strike would be 0 
    %for the GF function computation purpose of the origin would be the 
    %middle point of top edge
    X1=observation_points(:,1)*1e3-xorigin;
    Y1=observation_points(:,2)*1e3-yorigin;
    R1=[cosd(strike(k)) -sind(strike(k)); sind(strike(k)) cosd(strike(k))];
    Temp = R1*[X1';Y1'];
    X2=Temp(1,:)';Y2=Temp(2,:)';
    
    %To compute the GF, we need to refer to the middle point of the bootom
    %edge
    dbot =  z(k) + width(k)/2 .* sind(dip(k));
    shift = dbot/tand(dip(k)); 
    X=X2-shift; Y=Y2;
     [SSx1,SSy1,SSz]=green_function_okada(1,X,Y,nu,dip(k),dbot,length(k),...
        width(k),1,0);
     [DSx1,DSy1,DSz]=green_function_okada(1,X,Y,nu,dip(k),dbot,length(k),...
        width(k),2,0);
    
    %projection Based on the strike
    R2=[cosd(strike(k)) sind(strike(k)); -sind(strike(k)) cosd(strike(k))];
    Temp = R2*[SSx1';SSy1'];
    SSx=Temp(1,:);SSy=Temp(2,:);
    Temp = R2*[DSx1';DSy1'];
    DSx=Temp(1,:);DSy=Temp(2,:);
    
    %Build the Green function       
    G(1:n_surface_disp_dirs:end,k)= SSx;
    G(2:n_surface_disp_dirs:end,k)= SSy;
    G(3:n_surface_disp_dirs:end,k)= SSz;
    G(1:n_surface_disp_dirs:end,k+n_patches)= DSx;
    G(2:n_surface_disp_dirs:end,k+n_patches)= DSy;
    G(3:n_surface_disp_dirs:end,k+n_patches)= DSz;    
end



