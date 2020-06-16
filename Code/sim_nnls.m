function x = sim_nnls(A,d)
% FUCTION X = SIM_NNLS(A,D)
%
% This script is used to calculate the slip model by using LSQNONNEG method
% 
% by shwang @zhengzhou 2017-5-29
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
x = lsqnonneg(A,d);