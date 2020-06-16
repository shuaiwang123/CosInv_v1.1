function d = dist(x,y,x0,y0)
% x,y   : the centre position of triangles
% x0,y0 : the position of circle center
d = sqrt((x-x0).^2 + (y-y0).^2);
