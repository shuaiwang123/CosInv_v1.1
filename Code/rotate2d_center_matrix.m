function [x1,y1]=rotate2d_center_matrix(x,y,xc,yc,angle)
% this function is used to rotate a 2d point(x,y) of an
% angle(counterclockwise), the coordinates of rotation center is £¨xc,yc),
% the rotated point coordinates are the output (x1,y1)
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
[nc,nl] = size(x);

% the rotate matrix
R = [cosd(angle) -sind(angle);sind(angle) cosd(angle)];
%
x1 = zeros(size(x));
y1 = zeros(size(x));
%
for ii=1:nc
    for jj=1:nl
        rotX = R*[x(ii,jj)-xc;y(ii,jj)-yc];
        x1(ii,jj) = xc + rotX(1);
        y1(ii,jj) = yc + rotX(2);
    end
end
