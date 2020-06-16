function set_paths(code_dir)
%把code_dir目录下的文件夹和文件都添加到工作路径中
%
%...检查code_dir是否为目录...
if(~exist(code_dir,'dir'))
    error(['    [Directory of the Code Path',code_dir,'does not exist',...
        'Please check directory structure and spelling]']);
end
%...利用add_subdir_to_path()函数把code_dir都添加到工作路径中去...
add_subdir_to_path(code_dir);
disp(['    [Code directory has been added to path]']);

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%...Function add_subdir_to_patch...
function add_subdir_to_path(root_path)
addpath(root_path);            %把root_path添加到工作路径中去
files=dir(root_path);          %统计root_path下的文件夹和文件
for i=1:numel(files)
    if(files(i).isdir==1)
        if (~strcmp(files(i).name,'.') && ~strcmp(files(i).name,'..'))
            add_subdir_to_path([root_path,'/',files(i).name]);
        end
    end
end
%...

