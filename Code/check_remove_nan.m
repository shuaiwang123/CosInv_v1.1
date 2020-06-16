function temp=check_remove_nan(temp)
% this function is used to check whether the data line has NaN, if, then,
% remove it
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if ~isempty(temp(any(isnan(temp)'),:))
    disp('    [There is NaN in the dataline]');
    temp(any(isnan(temp)'),:) = [];
    temp = temp;
    disp('    [The NaN have been removed]');
end
