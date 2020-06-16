function [xy]= llh2localxy(llh,ll_origin)
% this function is used to convert lat-lon-height to local xy via a given
% origin, where x is the E, y is the North
%
% [lat;lon;height] to local coordinates xy=[x,y] based on the origin vector
% ll_origin=[lat_origin,long_origin], where x is the E, y is the North
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% the number of stations
[rows,nsta]=size(llh);
% change from deciaml degrees to decimal seconds
lat = 3600.0*llh(1,:);
lon = 3600.0*llh(2,:);

lat_origin = 3600*ll_origin(1);
%diff_lon   = 3600*ll_origin(2)*ones(size(lon)) - lon; % the original code
diff_lon   = 3600*ll_origin(2) - lon; % the code modified by ShuaiW

xy = zeros(nsta,2);
for i=1:nsta
    xy(i,:) = polyconic(lat(i),diff_lon(i),lat_origin);
end

% convert units from meter into kilometer and flip x-axis
xy(:,1) = -xy(:,1)/1000.0;  % the East
xy(:,2) =  xy(:,2)/1000.0;  % the North
