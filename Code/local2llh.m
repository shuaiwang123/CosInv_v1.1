function llh = local2llh(xy,origin)
% This function is used to convert the local coordinates to geographic
% coordinates.
%
% Input
%   XY : the local coordinates of points,the format is
%            1st         2nd            n
%            East        East           East
%            North       North          North
%    where,n is the number of points, unit gives in [km]
%   ORIGIN : the position of local coordinate's origin in geographic
%             coordinates, origin = [longitude latitude]
% Output
%    llh : the geographic coordinates of points, the format is 
%            1st         2nd            n
%            lon         lon            lon
%            lat         lat            lat
%   
% Note: this script is not a strict method to convert local coordiantes to
%       geographic coordinates, as this method do not take earth curvature
%       into cosideration, and the condition that 1 degree represents a
%       distance of 110km was assumed. We use this method for the sake of
%       plotting the dowansampled InSAR data in geographic coordinates with
%       the squares overlaid. There is a strict method to do this in the
%       ENIF document file.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
llh = repmat(origin(:),1,size(xy,2)) + xy./110;
