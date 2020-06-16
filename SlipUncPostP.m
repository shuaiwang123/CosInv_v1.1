% Slip uncertainties post-processing (SlipUncPostP)
%
% INPUT
%  MonteCarlo_filelist.txt its formats follows as
%  MonteCarlo_1_slip_model.mat
%  MonteCarlo_2_slip_model.mat
%  MonteCarlo_3_slip_model.mat
%           ...
%  CosInv_aketao_fault_trace_78_index6_linearramp.rec.mat
%  
%  Note that the last one is the optimal slip model
%
% REFERENCE
%  Qi, W et al.,(2011). Rupture of deep faults in the 2008 Wenchuan earthquake
%          and uplift of the Longmen Shan. Nature Geoscience, 4(9), 634-640.
%
% By shwang @whu 2017-03-18
% Modified Shuai WANG @POLYU 2018-03-16
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%
% close all;clear;clearvars;
% 
% filelist_name = 'MonteCarlo_filelist.txt';
% 
% % [ Load the slip mat file inverted from the perturbed noise datasets by 
% % [ the CosInvMonteCarlo_v1.0 ]
% filelist  = textread(filelist_name,'%s','headerlines',0);
% 
% % number of Monte Carlo datasets
% num_file  = numel(filelist);
% n_MC_datasets = num_file-1;
% 
% slip_tmp = {};
% for ii=1:num_file
%     filelist_tmp = char(filelist(ii));
%     load(filelist_tmp);
%     
%     slip_tmp{numel(slip_tmp)+1}  = slip;
% end
% 
% for ii=1:n_MC_datasets
%     slip_tmp{ii} = slip_tmp{ii} - slip_tmp{num_file};
% end
% 
% % [ Estimate the slip uncertainty ]
% slip_sigma = zeros(2*n_patches,1);
% 
% for ii=1:n_MC_datasets
%     slip_sigma = slip_sigma + slip_tmp{ii}.^2;
% end
% 
% slip = sqrt(slip_sigma/n_MC_datasets);
% 
% %clear slip_tmp slip_sigma 
% 
% save('MonteCarlo_slip_model_uncertainty.mat');


%% 
close all;clear;clearvars;

filelist_name = 'MonteCarlo_filelist.txt';

% [1. Load the slip mat file inverted from the perturbed noise datasets by the CosInvMonteCarlo_v1.0 ]
filelist  = textread(filelist_name,'%s','headerlines',0);

% number of Monte Carlo datasets
num_file  = numel(filelist);
n_MC_datasets = num_file;

for ii = 1
    filelist_ii = char(filelist(ii));
    load(filelist_ii);
end

% [2. Estimate the slip uncertainty ]
% matrix to store the Tslip for each MC 
Tslip_mc = zeros(n_patches,n_MC_datasets);


for mm=1:num_file
    filelist_ii = char(filelist(mm));
    load(filelist_ii);   
    Tslip_mc(:,mm) = Tslip(:);
end

slip_uncer_mc = std(Tslip_mc,[],2);

save MonteCarlo_slip_model_uncertainty.mat
        


