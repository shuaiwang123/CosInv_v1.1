function add_subdir_to_path(root_path)
% this function is used to add subdirectories to MATLAB path.
% adds root_path to the MATLAB path variable and then calls itself on 
% all subdirectories of ROOT_PATH.
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% add root_path to current MATLAB path
addpath(root_path);

% get list of all items contained in root_path 
files = dir(root_path);

% for each item of the list, check if it is a directory and if so call
% add_subdir_to_path on it
for i = 1:numel(files)
    if(files(i).isdir == 1) % if the file is a directory
        if(~ strcmp(files(i).name,'.') && ~ strcmp(files(i).name,'..')) % reject the directories '.' and '..'
            add_subdir_to_path([root_path,'/', files(i).name]); % call this script on the new directory
        end
    end
end