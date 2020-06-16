function plot_triangle_slip_distribution(fault_model,field,options)
% This function is used to plot slip distribution for triangle fault model
%
% 1.triangle fault model
%  fault_model = 
%  [x,y,z,strike,dip,area,vertices,strike_vect,updip_vect,normal_vect,rake] 
% 2.rectangle fault model
%  fault_model = 
%  [x,y,z,strike,dip,L,W,vertices,strike_vect,updip_vect,normal_vect,rake] 
%
% ### triangle plot options
% options = {...
%     'Perspective',ViewAngles,...
%     'FieldVector',slip_vector,vector_scale,vector_color,vector_width,...
%     'ColorMap',mycolormap,...
%     'AutoScale',1.4,...
%     'Figure',...
%     'Old',...
%     'ColorBar',...
%     'Shading','Flat',...
%     'ColorBarLabel','SLIP'};
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
x = fault_model(:,1);
y = fault_model(:,2);
z = fault_model(:,3);
if mean(mean(z)) < 0
    z = -z;
end

% [ The coordinates of three vetexes of a triangle elements ]
t1(:,1:3) = fault_model(:,7:9);
t2(:,1:3) = fault_model(:,10:12);
t3(:,1:3) = fault_model(:,13:15);

% The x,y,z coordinates of 1st vertex
xp = [t1(:,1) t2(:,1) t3(:,1)]';
% The x,y,z coordinates of 2nd vertex
yp = [t1(:,2) t2(:,2) t3(:,2)]';
% The x,y,z coordinates of 3nd vertex
zp = [t1(:,3) t2(:,3) t3(:,3)]';
up = [field(:) field(:) field(:)]';

% [ Default values ]
n_colorbar       = 0;
n_colorscale     = 0;
n_vector         = 0;
scale_vector     = 1;
n_colorbar_label = 0;
n_autoscale      = 0;
view_angles      = 2;
color_map        = jet;
symbol_size      = 5;
n_symbol         = 0;
shading_type     = 'Flat';
vect_lim         = 0.1;
gps_patchnumber  = 0;

% [ Load options ]
n_options = numel(options);
kk = 1;
while kk<=n_options
    switch options{kk}
        case 'Perspective'
            view_angles = options{kk+1};
            if ~isnumeric(view_angles)
                error('    [View angle not given properly]');
            end
            kk = kk + 1;
        case 'ColorMap'
            color_map = options{kk+1};
            kk = kk + 1;
        case 'ColorScale'
            color_scale_type = options{kk+1};
            if isnumeric(color_scale_type)
                n_colorscale = 1;
            else
                error('    [Color scale not given properly]');
            end
            kk = kk + 1;
        case 'FieldVector'
            field_vector = options{kk+1};
            if ~isnumeric(field_vector)
                error('    [Field vector not given properly]');
            end
            if isempty(find(size(field_vector) == 3))
                error('    [Field vector is not 3D]');
            end           
            n_vector = 1;
            kk = kk + 1;
            
            scale_vector = options{kk+1};
            if ~isnumeric(scale_vector)
                error('    [Vector scale not given properly]');
            end
            kk = kk + 1;
            
            vector_color = options{kk+1};
            if ~ischar(vector_color)
                error('    [Vector color not given properly]');
            end
            kk = kk + 1;
            
            vector_width = options{kk+1};
            if ~isnumeric(vector_width)
                error('    [Vector width not given properly]');
            end
            kk = kk + 1;
        case 'ColorBar'
            n_colorbar = 1;
        case 'AutoScale'
            n_autoscale = 1;
            auto_scale_factor = options{kk+1};
            if ~isnumeric(auto_scale_factor)
                error('    [Auto scale factor not given properly]');
            end
            kk = kk + 1;
        case 'PatchSymbol'
            n_symbol = 1;
            symbol_style = options{kk+1};
            if ~ischar(symbol_style)
                error('    [Symbol style not given properly]');
            end
            kk = kk + 1;
        case 'SymbolSize'
            symbol_size = options{kk+1};
            if ~isnumeric(symbol_size)
                error('    [Symbol size not given properly]');
            end
            kk = kk + 1;
        case 'Shading'
            shading_type = options{kk+1};
            if ~ischar(shading_type)
                error('    [Shading type not given properly]');
            end
            kk = kk + 1;
        case 'ColorBarLabel'
            n_colorbar_label = 1;
            color_bar_label = options{kk+1};
            if ~ischar(color_bar_label);
                error('    [Color bar not given properly]');
            end
            kk = kk + 1;
        case 'Gps_patchnumber'
            gps_patchnumber = options{kk+1};
            if ~isnumeric(gps_patchnumber)
                error('    [Gps_patchnumber not given properly]');
            end
            kk = kk + 1;
    end
    kk = kk + 1; % increment to next value
end

% [ Plotting ]
figure;
hold on;
colormap(color_map);
% 1. Plotting the center of fault patches
if gps_patchnumber == 1 
    if n_symbol == 1
        plot3(x,y,-z,symbol_style,'MarkeSize',symbol_size);
    end
    
% Note the fault number
    for i=1:size(x)
        text(x(i),y(i),-z(i)+1,[num2str(i)],'FontSize',10,'color','red');
    end
end

% 2. Plotting the slip distribution
patch(xp,yp,-zp,up);
% shading(shading_type);

% 3. Plotting the slip vector 
% Adjust the z direction of slip vector to meet plotting requirement 
if mean(mean(zp)) < 0
    zp = -zp;
    % Slip vector in geographic frame
    field_vector(:,3) = -field_vector(:,3);
end

% Only plotting the slip vector for slip amplitude > vect_lim*max(field)
icut = find(field > vect_lim*max(field));
if n_vector
    quiver3(x(icut),y(icut),-z(icut),field_vector(icut,1),field_vector(icut,2),...
        field_vector(icut,3),scale_vector,vector_color,'LineWidth',vector_width);
end

% 4. Plotting epicenter
scatter3(0,0,-15,300,'p','filled','r');

% 5. Set colorbar
if n_colorbar
    colorbar('location','EastOutside');
    caxis([0 0.99*max(field)]);
    if n_colorbar_label
        colorbar('FontSize',20),set(get(colorbar,'ylabel'),'string',color_bar_label);
    end
end
view(view_angles);

% 6. Set the axis range
if n_autoscale
    fx_figure = auto_scale_factor;
    fy_figure = fx_figure;
    xmin = fx_figure*min(fault_model(:,1));
    xmax = fx_figure*max(fault_model(:,1));
    ymin = fy_figure*min(fault_model(:,2));
    ymax = fy_figure*max(fault_model(:,2));
    zmin = -max(fault_model(:,3))-5;
    zmax = 0; 
    axis([ xmin xmax ymin ymax zmin zmax]);
end 

% 7. Plotting settings
% axis equal;
axis on;grid on;box on;

% Set tick 
set(gca,'FontSize',20); 
set(gca,'xminortick','off');
set(gca,'yminortick','off');
set(gca,'zminortick','off');
set(gca,'ticklength',[0.01 0.01]);
set(gca,'tickdir','out');

% Set tick label
% % --- Synthetic Tests ---
% set(gca,'XLim',[-95 95]);
% set(gca,'XTick',-80:20:80);
% set(gca,'XTickLabel',{'-80','-60','-40','-20','0','20','40','60','80'});
% set(gca,'YLim',[0 160]);
% set(gca,'YTick',0:20:140);
% set(gca,'YTickLabel',{'0','20','40','60','80','100','120','140','160'});
% set(gca,'ZLim',[-35 0]);
% set(gca,'ZTick',-30:10:0);
% set(gca,'ZTickLabel',{'-30','-20','-10','0'});

% --- Nepal Mw7.8 ---
% set(gca,'XLim',[-55 195]);
% set(gca,'XTick',-50:25:175);
% set(gca,'XTickLabel',{'-50','-25','0','25','50','75','100','125','150','175'});
% set(gca,'YLim',[-140 125]);
% set(gca,'YTick',-130:30:110);
% set(gca,'YTickLabel',{'-130','-100','-70','-40','-10','20','50','80','110'});
% set(gca,'ZLim',[-35 0]);
% set(gca,'ZTick',-30:10:0);
% set(gca,'ZTickLabel',{'30','20','10','0'});

% --- Aketao Mw6.7 ---
% set(gca,'XLim',[-30 70]);
% set(gca,'XTick',-30:25:70);
% set(gca,'XTickLabel',{'-30','-5','20','45','70'});
% set(gca,'YLim',[-50 20]);
% set(gca,'YTick',-50:20:20);
% set(gca,'YTickLabel',{'-50','-30','-10','10'});
% set(gca,'ZLim',[-35 0]);
% set(gca,'ZTick',-30:10:0);
% set(gca,'ZTickLabel',{'30','20','10','0'});

% Set label
xlabel('East (km)','fontsize',20,'FontName','Arial','FontWeight','normal','rotation',14);
ylabel('North (km)','fontsize',20,'FontName','Arial','FontWeight','normal','rotation',-14);
zlabel('Depth (km)','fontsize',20,'FontName','Arial','FontWeight','normal');
title('Triangle fault patches');

view(view_angles);
xlabel('East (km)');
ylabel('North (km)');
zlabel('Depth (km)');









