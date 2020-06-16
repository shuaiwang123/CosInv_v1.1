function slip_vector = build_slip_vectors(U_str,U_updip,fault_model)
% This function is used to take U_str and U_updip on the local fault axes,
% and tranform it as a unit vector in geographical(east,north,up) reference
% frame based on fault elements defined in FAULT_MODEL
%
% 1.triangle fault model
%  fault_model = 
%  [x,y,z,strike,dip,area,vertices,strike_vect,updip_vect,normal_vect,rake] 
% 2.rectangle fault model
%  fault_model = 
%  [x,y,z,strike,dip,L,W,vertices,strike_vect,updip_vect,normal_vect,rake] 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% The fixed number of parameters for rectangle fault model
n_param_rect = 29;
% The fixed number of parameters for point/triangle fault model
n_param_trg  = 25;

if size(fault_model,2) == n_param_rect
    ind_strike_vect = 20:22;
    ind_updip_vect  = 23:25;
elseif size(fault_model,2) == n_param_trg
    ind_strike_vect = 16:18;
    ind_updip_vect  = 19:21;
else
    error('    [The fault model does not correspond neither to a rectangle model nor triangle model]');
end

n_patches = numel(U_str);
% Initialization of output vectors
slip_vector = zeros(n_patches,3);

% [ Convert slip vector on fault to geographical frame ]
for ii=1:n_patches
    Vstrike = fault_model(ii,ind_strike_vect);
    Vupdip  = fault_model(ii,ind_updip_vect);
    
    slip_vector(ii,:) = U_str(ii)*Vstrike + U_updip(ii)*Vupdip;
    slip_vector(ii,:) = slip_vector(ii,:)/norm(slip_vector(ii,:),'fro');
end

    
    
    
    
    
    
    


    