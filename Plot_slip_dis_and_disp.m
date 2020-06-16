% This scripts is used to plot RECTANGLE slip distribution.
%
% Note that the unit of the inverted slip and synthetic and simulated
% displacement used in the plot is in [m] by default. So if the observed
% data used to invert the slip is not in [m], before plotting you should
% transform them to [m]
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
disp('...... Plotting ......');
% =========================================================================
% ====                        Plotting Options                         ====
% =========================================================================
% Define the view angles using median strike and dip
AZ      = 3;%133;
EL      = 18;%44;
ViewAngles = [AZ,EL];

% Defining the color map
mycolormap = colormap(flipud(hot));
% mycolormap = colormap(cmap);
% mycolormap = colormap(flipud(gray));
% cmap       = load('Code/Plotting/precip2_17lev2matlabColorbar.cpt');
mycolormap = colormap(mycolormap);

% [ Build model data ]
% synthetic GPS 
syngps  = dmodel(d_index{1});
syngpse = zeros(gps.nobs,1);
syngpsn = zeros(gps.nobs,1);
syngpsu = zeros(gps.nobs,1);
indexe  = 1;
indexn  = 2;
indexu  = 3;
for ii = 1:gps.nobs
    syngpse(ii) = syngps(indexe);
    syngpsn(ii) = syngps(indexn);
    syngpsu(ii) = syngps(indexu);
    indexe      = indexe + 3;
    indexn      = indexn + 3;
    indexu      = indexu + 3;
end

% [ Build the slip vector ]
% strike-slip
slip_x = slip(1:n_patches,1)*UNIT2m;
% dip-slip
slip_y = slip(n_patches+1:end,1)*UNIT2m;
% total slip
Tslip  = sqrt(slip_x.^2 + slip_y.^2);
% build slip vector
slip_vector = build_slip_vectors(slip_x,slip_y,fault_model);
% options to plot the slip vector
vector_scale = 0.5;
vector_color = 'w';
vector_width = 1;
% options to plot GPS and patch number 
gps_patchnumber = 0;%input('    [Please choosing whether plotting the station name and disp vectors of GPS and the patch number: 0. NO  1.YES]');

% =========================================================================
% ====   Plotting slip distribution and observed data and model data   ====
% =========================================================================
% [ Varous options for plot_rectangle_slip_distribution_and_station_vector ]
plot_options = {...
    'Perspective',ViewAngles,...
    'FieldVector',slip_vector,vector_scale,vector_color,vector_width,...
    'ColorMap',mycolormap,...
    'AutoScale',1,...
    'Gps_patchnumber',gps_patchnumber,...
    'Figure',...
    'Old',...
    'ColorBar',...
    'Shading','Flat',...
    'ColorBarLabel',['Slip (m)']};

% [ Plotting slip distribution ]
if size(fault_model,2) == 29
    plot_rec_slip_dis(fault_model,Tslip,plot_options);
else
    plot_trg_slip_dis(fault_model,Tslip,plot_options);
end

% [ Plotting observed and model data ]
disp('    [Plotting observed data and model data]');
% plotting GPS data
figure;
subplot(1,2,1);
box on;axis square;
quiver(gps.x,gps.y,gps.disp(:,1)*UNIT2cm,gps.disp(:,2)*UNIT2cm,'k','linewidth',6);hold on;
quiver(gps.x,gps.y,syngpse*UNIT2cm,syngpsn*UNIT2cm,'r');
xlabel('East(km)');ylabel('North(km)');title('Horizontal GPS');
subplot(1,2,2);
box on;axis square;
quiver(gps.x,gps.y,zeros(gps.nobs,1),gps.disp(:,3)*UNIT2cm,'k','linewidth',6);hold on;
quiver(gps.x,gps.y,zeros(gps.nobs,1),syngpsu*UNIT2cm,'r');
xlabel('East(km)');ylabel('North(km)');title('Vertical GPS');

% plotting InSAR data
cmap = load('Code/Plotting/bycr2matlabColorbar.cpt');
for ii=2:n_datasets
    figure;
    subplot(1,3,1);
    scatter(x(position_index{ii}),y(position_index{ii}),10,d(d_index{ii})*UNIT2cm,'filled');
    box on;axis on;axis image;colormap(cmap);set(get(colorbar,'ylabel'),'string','LOS(cm)');
    xlabel('East(km)');ylabel('North(km)');title('Observed InSAR');
%     if ii == 2,
%         caxis([-20 100]);
%     elseif ii == 3,
%         caxis([-70 75]);
%     end
%     caxis([-10 10]);
%     xlim([-20 40]);ylim([-40 30]);
    
    subplot(1,3,2);
    scatter(x(position_index{ii}),y(position_index{ii}),10,dmodel(d_index{ii})*UNIT2cm,'filled');
    box on;axis on;axis image;colormap(cmap);set(get(colorbar,'ylabel'),'string','LOS(cm)');
    xlabel('East(km)');ylabel('North(km)');title('Simulated InSAR');
%     if ii == 2,
%         %caxis([-20 100]);
%     elseif ii == 3,
%         %caxis([-70 75]);
%     end
%     caxis([-10 10]);
%     caxis([-10 10]);
%     xlim([-20 40]);ylim([-40 30]);

    subplot(1,3,3);
    scatter(x(position_index{ii}),y(position_index{ii}),10,d(d_index{ii})*UNIT2cm - dmodel(d_index{ii})*UNIT2cm,'filled');
    box on;axis on;axis image;colormap(cmap);set(get(colorbar,'ylabel'),'string','LOS(cm)');
    xlabel('East(km)');ylabel('North(km)');title('Residual InSAR');
%     if ii == 2,
%         %caxis([-20 100]);
%     elseif ii == 3,
%         %caxis([-70 75]);
%     end
%     caxis([-10 10]);
%     caxis([-10 10]);
%     xlim([-20 40]);ylim([-40 30]);
end









