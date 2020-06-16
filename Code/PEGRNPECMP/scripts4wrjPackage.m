% This script is used to prepare the RECTANGULAR SUBFAULTS parameters
% required by GFZ Prof. Wang package. 
%
% By Shuai WANG @POLY 2017-12-31, the last day of 2017, new year is coming
%
%==========================================================================

%% Prepare pegrn20160520-input.dat
pegrninp = 'pegrn20160520-input.dat';
npatches = 340;
np_st    = 34;
np_di    = 10;
patsize  = 0.8;

% Remove the old pegrninp file and building a new one
if exist(pegrninp,'file'),delete(pegrninp);end

% determine the pos_s and pos_d
pos_s = zeros(npatches,1);
pos_d = zeros(npatches,1);
for ii = 1:np_st
    for jj = 1:np_di
        pos_s(jj + (ii-1)*np_di) = patsize/2 + (ii-1)*patsize;
        pos_d(jj + (ii-1)*np_di) = patsize/2 + (jj-1)*patsize;
    end
end    
     
if ~exist(pegrninp,'file')
    fid = fopen(pegrninp,'w');
    
    fprintf(fid,'# n_faults\n');
    fprintf(fid,'#-------------------------------------------------------------------------------\n');   
    fprintf(fid,'%s\n','1');
    fprintf(fid,'#-------------------------------------------------------------------------------\n');
    fprintf(fid,'# n   O_lat   O_lon    O_depth length  width strike dip   np_st np_di start_time\n');
    fprintf(fid,'# [-] [deg]   [deg]    [km]    [km]     [km] [deg]  [deg] [-]   [-]   [day]\n');
    fprintf(fid,'#     pos_s   pos_d    slp_stk slp_dip open\n');
    fprintf(fid,'#     [km]    [km]     [m]     [m]     [m]\n');
    fprintf(fid,'#-------------------------------------------------------------------------------\n');
    fprintf(fid,'%4i%10.4f%10.4f%10.4f%10.4f%10.4f%15.4f%10.4f%4i%4i%5.2f\n',   1,    -25.6959,   129.9465,     0,    27,    8,    312,     22,   34,   10,   0);
  
    for ii = 1:npatches       
        fprintf(fid,'%14.4f%10.4f%10.4f%10.4f%10.4f\n',    pos_s(ii),    pos_d(ii),         slip_x(ii),         -slip_y(ii),         0);
    end
end
