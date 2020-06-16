% function plotCFS(XGRID,YGRID,coulomb,titlename)
function plotShadeMap(XGRID,YGRID,coulomb,mask_file_name,titlename)
% PLOTCFS
%   This function can plot SHEAR/NORMAL STRESS and CFS
%
% INPUT
%   XGRID,YGRID    the CFS grid points
%   coulomb        CFS[Bar]
%   mask_file_name is usd to extract pixels falled into a poly
%   titlename      the figure you want to plot, can be SHEAR/NORMAL/CFS
%
% Modifed from plotCFS
% by shwang @whu 2017-05-22
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% set(gca,'FontSize',20); 
% hold on;box on;
[X,Y,Z] = griddata(XGRID,YGRID,coulomb,linspace(min(XGRID),max(XGRID),250)',linspace(min(YGRID),max(YGRID),250),'v4'); % 'v4'/'cubic'/'nearest'
% pcolor(X,Y,Z);shading interp;axis image;

X = X(:);
Y = Y(:);
Z = Z(:);

mask_file = load(mask_file_name);
mask_lon  = mask_file(:,1);
mask_lat  = mask_file(:,2);

% disp('    [Extract Stations]');
in    = inpolygon(X,Y,mask_lon,mask_lat);
index = find(in ~= 0);

X     = X(index);
Y     = Y(index);
Z     = Z(index);


scatter(X,Y,5,Z,'filled');

%axis([84.5 86 27.5 29]);
% caxis([-1,1]);
% colorbar('FontSize',20);set(get(colorbar,'ylabel'),'string','Bar');
colorbar;
% xlabel('East (km)','fontsize',20,'FontName','Arial','FontWeight','normal','rotation',0);
% ylabel('North (km)','fontsize',20,'FontName','Arial','FontWeight','normal','rotation',90);
title(titlename);
% hold off;