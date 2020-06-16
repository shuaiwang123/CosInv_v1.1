function plot_rec_slip_dis(fault_model,field,options)
% This function is used to plot slip distribution for rectangle fault model
%
% 1.triangle fault model
%  fault_model = 
%  [x,y,z,strike,dip,area,vertices,strike_vect,updip_vect,normal_vect,rake] 
% 2.rectangle fault model
%  fault_model = 
%  [x,y,z,strike,dip,L,W,vertices,strike_vect,updip_vect,normal_vect,rake] 
%
% ### rectangle plot options
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
disp('    [Plotting rectangle source slip distribution]');
x = fault_model(:,1);               % E-W. Unit[km]
y = fault_model(:,2);               % N-S. Unit[km]
z = fault_model(:,3);               % Down is positive. Unit[km]
strike = fault_model(:,4);
dip    = fault_model(:,5);
dl     = fault_model(:,6);
dw     = fault_model(:,7);
if mean(mean(z)) < 0
    z = -z;
end

% [ Default values ]
n_colorbar       = 0;
n_colorscale     = 0;
n_vector         = 0;
scale_vector     = 0;%1;
n_colorbar_label = 0;
n_autoscale      = 0;
view_angles      = 2;
color_map        = jet;
symbol_size      = 5;
n_symbol         = 0;
shading_type     = 'Flat';
vect_lim         = 0.1;
gps_patchnumber  = 1;

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
index_vertices = 8:19;
vertices = fault_model(:,index_vertices);
xp = [vertices(:,1) vertices(:,4) vertices(:,7) vertices(:,10)]';
yp = [vertices(:,2) vertices(:,5) vertices(:,8) vertices(:,11)]';
zp = [vertices(:,3) vertices(:,6) vertices(:,9) vertices(:,12)]';

slip_temp = [field(:) field(:) field(:) field(:)]';
patch(xp,yp,-zp,slip_temp);
% shading(shading_type);

% 3. Plotting epicenter
% === 2016/11/25 Mw6.7 Aketao earthquake ===
% scatter3(0,-5,-9.703,300,'p','filled','w');
% === 2009/01/03 Mw7.6 Irian earthquake ===
% % Poiata 2010 EPS
% scatter3(0,0,-10,300,'p','filled','r');
% scatter3(62.33,-22.11,-12,300,'p','filled','r');
% % scatter3(0,0,0,300,'p','filled','m');     % checkerboard
% % scatter3(62.33,-22.11,0,300,'p','filled','m');

% GCMT
% scatter3(6.5657,14.3737,-15.2,300,'p','filled','r');
% scatter3(77.9152,-7.7397,-18.2,300,'p','filled','r');
% === 2017/05/29 Mw 6.6 Indonesia earthquake ===
scatter3(0,0,0,300,'p','filled','r');

if gps_patchnumber == 0
%     origin      = [39.27 74.04];
%     aftershocks = load('E:\Coseismic\CosInv_v1.0\Irian Coseismic Disp\aftershocks_proj_Fang.txt');
% %     lon_after   = aftershocks(:,1);
% %     lat_after   = aftershocks(:,2);
% %     dep_after   = aftershocks(:,3);
% %     llh_after   = [lat_after';lon_after'];
% %     xy_after    = llh2localxy(llh_after,origin);
% %     x_after     = xy_after(:,1);
% %     y_after     = xy_after(:,2);
% %     z_after     = dep_after;
%     x_after   = aftershocks(:,1);
%     y_after   = aftershocks(:,2);
%     z_after   = aftershocks(:,3);
%     scatter3(x_after,y_after,z_after,80,'o','k');
end

% 4. Plotting the slip vector 
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

% 5. Set colorbar
if n_colorbar
    colorbar('location','EastOutside');
    caxis([0 0.99*max(field)]);
    if n_colorbar_label
        colorbar('FontSize',20),set(get(colorbar,'ylabel'),'string',color_bar_label);
    end
end

% 6. Set the axis range
if n_autoscale
%     fx_figure = auto_scale_factor;
%     fy_figure = fx_figure;
%     xmin = fx_figure*min(fault_model(:,1));
%     xmax = fx_figure*max(fault_model(:,1));
%     ymin = fy_figure*min(fault_model(:,2));
%     ymax = fy_figure*max(fault_model(:,2));
%     axis([ xmin xmax ymin ymax]);
end 

% 7. Plotting settings
% axis equal;
axis on;box on;% grid on;
caxis([0 1]);axis image;

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
% set(gca,'XLim',[-40 60]);
% set(gca,'XTick',-40:20:60);
% set(gca,'XTickLabel',{'-40','-20','0','20','40','60'});
% set(gca,'YLim',[-40 10]);
% set(gca,'YTick',-40:20:10);
% set(gca,'YTickLabel',{'-40','-20','0'});
% set(gca,'ZLim',[-35 0]);
% set(gca,'ZTick',-30:10:0);
% set(gca,'ZTickLabel',{'30','20','10','0'});

% --- Irian Mw7.6 ---
% set(gca,'XLim',[-120 150]);
% set(gca,'XTick',-80:40:120);
% set(gca,'XTickLabel',{'-80','-40','0','40','80','120'});
% set(gca,'YLim',[-70 70]);
% set(gca,'YTick',-70:35:70);
% set(gca,'YTickLabel',{'-70','-35','0','35','70'});
% set(gca,'ZLim',[-50 0]);
% set(gca,'ZTick',-50:25:0);
% set(gca,'ZTickLabel',{'50','25','0'});

% % --- Indonesia Mw 6.6 ---
% set(gca,'XLim',[-25 25]);
% set(gca,'XTick',-25:25:25);
% set(gca,'XTickLabel',{'-25','0','25'});
% set(gca,'YLim',[-25 15]);
% set(gca,'YTick',-20:10:10);
% set(gca,'YTickLabel',{'-20','-10','0','10'});
% set(gca,'ZLim',[-30 0]);
% set(gca,'ZTick',-30:15:0);
% set(gca,'ZTickLabel',{'30','15','0'});


% --- Papua Mw 7.5 ---
% set(gca,'XLim',[-100 200]);
% set(gca,'XTick',-100:50:200);
% set(gca,'XTickLabel',{'-100','-50','0','50','100','150','200'});
% set(gca,'YLim',[-100 100]);
% set(gca,'YTick',-100:50:100);
% set(gca,'YTickLabel',{'-100','-50','0','50','100'});
% set(gca,'ZLim',[-40 0]);
% set(gca,'ZTick',-40:20:0);
% set(gca,'ZTickLabel',{'40','20','0'});

% Set label
xlabel('East (km)','fontsize',20,'FontName','Arial','FontWeight','normal','rotation',0);
ylabel('North (km)','fontsize',20,'FontName','Arial','FontWeight','normal','rotation',45);
zlabel('Depth (km)','fontsize',20,'FontName','Arial','FontWeight','normal');
title('Rectangle fault patches');

view(view_angles);



    
    


            












