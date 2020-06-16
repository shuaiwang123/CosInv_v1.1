% This script is used to load fault model including rectangle disloacation
% model, triangle dislocation model and iregular dislocation model.
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
disp('...... Load Fault Model ......');

% ============ Synthetic Fault Model ============
% fault_model_file = 'model_iregular.trg';
% fault_model_file = 'model_curved.trg';
fault_model_file = 'model_curved.rec';
% fault_model_file = 'model.trg';
% fault_model_file = 'model.rec';

fault_model = load(fault_model_file);
if size(fault_model,2) == 29
    disp('    [RDEs model have been loaded]');
    xc          = fault_model(:,1);
    yc          = fault_model(:,2);
    zc          = fault_model(:,3);
    strike      = fault_model(:,4);
    dip         = fault_model(:,5);
    leng        = fault_model(:,6);
    width       = fault_model(:,7);
    vertices    = fault_model(:,8:19);
    strike_vect = fault_model(:,20:22);
    updip_vect  = fault_model(:,23:25);
    normal_vect = fault_model(:,26:28);
    rake        = fault_model(:,29);   
    
    fault_model = [xc,yc,zc,strike,dip,leng,width,vertices,strike_vect,updip_vect,normal_vect,rake];
    
elseif size(fault_model,2) == 25
    disp('    [TDEs model have been loaded]');
    xc           = fault_model(:,1);
    yc           = fault_model(:,2);
    zc           = fault_model(:,3);
    strike       = fault_model(:,4);
    dip          = fault_model(:,5);
    area         = fault_model(:,6);
    vertices     = fault_model(:,7:15);
    strike_vect  = fault_model(:,16:18);
    updip_vect   = fault_model(:,19:21);
    normal_vect  = fault_model(:,22:24);
    rake         = fault_model(:,25);
    
    fault_model  = [xc,yc,zc,strike,dip,area,vertices,strike_vect,updip_vect,normal_vect,rake];
end

    
    


