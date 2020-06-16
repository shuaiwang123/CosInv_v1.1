% This script is used to save the CosInv_v1.0 inversion results.
%
% The format of output file is
%  XC[km]  YC[km]  LON[Degree] LAT[Degree]  STR_SLIP[m]  REV_SLIP[m]  SLIP_VECTOR  TSLIP[m]
%
%  XC,YC          local coordinates of fault center
%  LON,LAT        geographic coordinates of fault center
%  STR_SLIP       strike slip component
%  REV_SLIP       reverse slip component
%  TSLIP          total slip of fault patches
%  SLIP_VECTOR    slip vector in geographic frame
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
disp('...... Save Inverted Results ......');
% [ Convert the local coordinates of fault to geographic coordinates ]
xy_fault  = [fault_model(:,1)';fault_model(:,2)'];
llh_fault = local2llh(xy_fault,origin);

% [ Ouput the inversion results]
% cd Output_files;
% header1 = ' xc      yc        zc        lon      lat      str_slip      thrus_slip      total_slip      slip_vect_e      slip_vect_n      slip_vect_n';
% header2 =' [km]    [km]       [km]      [Deg]     [Deg]      [m]            [m]             [m]             -                 -                 -';             
% output_inver_result = [fault_model(:,1),fault_model(:,2),fault_model(:,3),llh_fault(1,:)',llh_fault(2,:)',slip_x,slip_y,Tslip,slip_vector];
header1 = ' lon        lat       zc      strike_slip      thrust_slip     ';
header2 =' [Deg]      [Deg]     [km]       [m]                [m]         ';             
output_inver_result = [llh_fault(1,:)',llh_fault(2,:)',fault_model(:,3),slip_x,slip_y];

output_inver_result = check_remove_nan(output_inver_result);

dlmwrite([fault_model_file,'.slip'],header1,'delimiter','');
dlmwrite([fault_model_file,'.slip'],header2,'-append','delimiter','');
dlmwrite([fault_model_file,'.slip'],output_inver_result,'-append','delimiter','\t','precision','%.4f');

% [ Output the simulated displacement ]
% % save simulated gps
% header1        = '% lon     lat    ve   vn    up';
% header2        = '% [Deg]  [Deg]  [the same unit as the input]';
% output_sim_gps = [gps.lon,gps.lat,syngpse,syngpsn,syngpsu];
% dlmwrite(['sim_gps_',fault_model_file],header1,'delimiter','');
% dlmwrite(['sim_gps_',fault_model_file],header2,'-append','delimiter','');
% dlmwrite(['sim_gps_',fault_model_file],output_sim_gps,'-append','delimiter','\t','precision','%.4f');
% 
% % save simulated insar
% for ii=2:n_datasets 
%     header1                = '% lon     lat               los';
%     header2                = '% [Deg]  [Deg]   [the same unit as the input]';
%     output_sim_insar(ii-1) = [lon(d_index{ii}),lat(d_index{ii}),dmodel(d_index{ii})];
%     dlmwrite(['sim_insar',num2str(ii-1),'_',fault_model_file],header1,'delimiter','');
%     dlmwrite(['sim_insar',num2str(ii-1),'_',fault_model_file],header2,'-append','delimiter','');
%     dlmwrite(['sim_insar',num2str(ii-1),'_',fault_model_file],output_sim_insar(ii-1),'-append','delimiter','\t','precision','%.4f');
% end

% [ Output the inversion results that required by CFS_v1.0 ]
if size(fault_model,2) == 29
    header1 = '% vertices(1~12)   strike    dip      str_slip    thrus_slip     tensile_slip';
    header2 = '%    [km]           [Deg]    [Deg]     [m]            [m]             [m]';
    output4CFS = [fault_model(:,8:19),fault_model(:,4),fault_model(:,5),slip_x,slip_y,zeros(numel(slip_x),1)];
    output4CFS = check_remove_nan(output4CFS);
    dlmwrite([fault_model_file,'.inp'],header1,'delimiter','');
    dlmwrite([fault_model_file,'.inp'],header2,'-append','delimiter','');
    dlmwrite([fault_model_file,'.inp'],output4CFS,'-append','delimiter','\t','precision','%.4f');
elseif size(fault_model,2) == 25
    header1 = '% xc    yc    zc    vertices(4~12)   str_slip    thrus_slip     tensile_slip';
    header2 = '%[km]  [km]  [km]     [km]            [m]            [m]             [m]';
    output4CFS = [fault_model(:,1:3),fault_model(:,7:15),slip_x,slip_y,zeros(numel(slip_x),1)];
    output4CFS = check_remove_nan(output4CFS);
    dlmwrite([fault_model_file,'.inp'],header1,'delimiter','');
    dlmwrite([fault_model_file,'.inp'],header2,'-append','delimiter','');
    dlmwrite([fault_model_file,'.inp'],output4CFS,'-append','delimiter','\t','precision','%.4f');
end

% [ Output the inversion results that required by SimDisp_v1.0 ]
if size(fault_model,2) == 29 
    header1 = '%xc    yc    zc  strike   dip   leng  width  slip_str  slip_dip';
    header2 = '%[km]  [km] [km]  [Deg]  [Deg]  [km]  [km]     [m]        [m]';
    output4SimDisp = [fault_model(:,1:7),slip_x,slip_y];
    output4SimDisp = check_remove_nan(output4SimDisp);
    dlmwrite([fault_model_file,'4SimDisp.slip'],header1,'delimiter','');
    dlmwrite([fault_model_file,'4SimDisp.slip'],header2,'-append','delimiter','');
    dlmwrite([fault_model_file,'4SimDisp.slip'],output4SimDisp,'-append','delimiter','\t','precision','%.4f');
elseif size(fault_model,2) == 25
    header1 = '%xc    yc    zc  strike   dip   area  vertices slip_str  slip_dip';
    header2 = '%[km]  [km] [km]  [Deg]  [Deg]  [km2]    -        [m]        [m]';
    output4SimDisp = [fault_model(:,1:15),slip_x,slip_y];
    output4SimDisp = check_remove_nan(output4SimDisp);
    dlmwrite([fault_model_file,'4SimDisp.slip'],header1,'delimiter','');
    dlmwrite([fault_model_file,'4SimDisp.slip'],header2,'-append','delimiter','');
    dlmwrite([fault_model_file,'4SimDisp.slip'],output4SimDisp,'-append','delimiter','\t','precision','%.4f');
end

save(['CosInv_',fault_model_file,'.mat']);
%cd ..;
