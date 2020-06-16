% This scripts is used to invert distributed slip at depth
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
disp('...... Invert for Distributed Slip at Depth ......');


% [ Transform the input unit to m ]
UNITinput = 'cm';
% for plotting inverted slip and simulated displacement
UNIT2m    = 0;    
% for plotting Lcurve
UNIT2cm   = 0;
switch UNITinput
    case 'mm', UNIT2m = 10^-3; UNIT2cm = 10^-1;
    case 'cm', UNIT2m = 10^-2; UNIT2cm = 10^0;
    case 'dm', UNIT2m = 10^-1; UNIT2cm = 10^1;
    case 'm',  UNIT2m = 10^0;  UNIT2cm = 10^2;
end

% [ LAP weight ]
% ============ Synthetic Test ============
% lap_weight = 10.^([linspace(-2,0.6,10)]);% model_iregualr.trg
lap_weight = 10.^([linspace(-0.5,0.5,10)]);  % model_curved.trg/model_curved.rec/model.trg
% lap_weight = 10.^([linspace(-2,1.5,10)]);
% lap_weight = 10.^([linspace(-0.5,1,10)]);

% [ LSQLIN Options ]
% ============ Synthetic Test ============
% A*mil < B
A       = blkdiag(0*eye(n_patches),  -eye(n_patches),  zeros(nramp)); % Forces invser slip and allows right and left leteral slip
B       = zeros(2*n_patches+nramp,  1); 
% lb < mil < ub
lb      = [ -0.1*ones(n_patches,1);  0*ones(n_patches,1);   -9999*ones(nramp,1)]; 
ub      = [  0.1*ones(n_patches,1);  600*ones(n_patches,1);    9999*ones(nramp,1)];

% % Initial value 
% s0      = 0.5*(lb + ub);
% options = optimset('LargeScale','off','MaxIter',4000);
% % options = optimoptions('LargeScale','off','MaxIter',4000);

% [ Construct Design matrix and Observed vector ]
ch      = chol(covd);
Icovd   = inv(ch');
G       = Icovd * [G_proj rampg];   % Weighted GREENs function
D       = [Icovd*d; zeros(2*n_patches,1)]; % Weighted observed data

% [ Inversion ]
n_lap        = length(lap_weight);
missfit      = zeros(n_lap,1);
roughness    = zeros(n_lap,1);
rmsgps_insar = cell(n_datasets,1);

method_flag = input('    [Please choose the method to invert slip: 1. Classical LSQ  2.BVLS  3. LSQLIN]');
for i=1:n_lap
    lambda  = lap_weight(i);
    Gsmooth = [G; lambda*Lap zeros(2*n_patches,nramp)]; % Design matrix for inversion
    
    % Initial value
    %s0      = 0.5*(lb + ub);
    options = optimset('LargeScale','off','MaxIter',2000);
    % options = optimoptions('LargeScale','off','MaxIter',4000);

    if method_flag == 1 ,disp('      ... Using Classical LSQ to invert for slip ...');
        s = (Gsmooth'*Gsmooth)^-1 * Gsmooth' * D;
    elseif method_flag ==2, disp('      ... Using BVLS to invert for slip ...');
        s0     = sim_nnls(Gsmooth, D);
        opt.x0 = s0;
        s      = sim_bvls(Gsmooth, D, lb, ub, opt);
    else  disp('      ... Using LSQLIN to invert for slip ...');
        Gsmooth = double(Gsmooth);
        D       = double(D);
        A       = double(A);
        B       = double(B);
        lb      = double(lb);
        ub      = double(ub);
        s0      = sim_nnls(Gsmooth, D);
        s0      = double(s0);
        [s,resnorm,residual] = lsqlin(Gsmooth, D, A, B, [], [], lb, ub, s0, options);
    end
    
    dmodel    = ch' * G * s;
    dres      = d - dmodel; % Residual between model data and observed data
    %rms       = (dres(:)'*dres(:)/numel(d))^0.5;
    missfit(i)= dres(:)'*dres(:);
    for j=1:n_datasets
        rmsgps_insar{j} = (dres(d_index{j})'*dres(d_index{j})/numel(d_index{j}))^0.5;
    end
    % 
    % roughness(i) = sum(abs([Lap zeros(2*n_patches,nramp)]*s))/(2*n_patches);
    roughness(i) = sum(abs([Lap zeros(2*n_patches,nramp)]*s))/n_patches;
end

% [ Choosing the Best Laplacian weight using L_curve ]
lambda  = plot_lcurve(missfit,roughness,lap_weight,UNIT2cm);
Gsmooth = [G; lambda*Lap zeros(2*n_patches,nramp)];
Resol   = inv(Gsmooth'*Gsmooth)*G'*G; % Resolution matrix of model

% [ Inversion using the best Laplacian weight ]
if method_flag == 1 
    s = (Gsmooth'*Gsmooth)^-1 * Gsmooth' * D;
else  
    [s,resnorm,residual] = lsqlin(Gsmooth, D, A, B, [], [], lb, ub, s0, options);
end
% dmodel = G*m + L*s (ramp)
dmodel = ch' * G * s;
slip   = s(1:2*n_patches);
    
    
    
    
    