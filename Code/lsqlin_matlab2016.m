function [X,resnorm,residual,exitflag,output,lambda] = lsqlin(C,d,A,b,Aeq,beq,lb,ub,X0,options,varargin)
%LSQLIN Constrained linear least squares.
%   X = LSQLIN(C,d,A,b) attempts to solve the least-squares problem
%
%           min  0.5*(NORM(C*x-d)).^2       subject to    A*x <= b
%            x
%
%   where C is m-by-n.
%
%   X = LSQLIN(C,d,A,b,Aeq,beq) solves the least-squares
%   (with equality constraints) problem:
%
%           min  0.5*(NORM(C*x-d)).^2    subject to 
%            x                               A*x <= b and Aeq*x = beq
%
%   X = LSQLIN(C,d,A,b,Aeq,beq,LB,UB) defines a set of lower and upper
%   bounds on the design variables, X, so that the solution
%   is in the range LB <= X <= UB. Use empty matrices for 
%   LB and UB if no bounds exist. Set LB(i) = -Inf if X(i) is unbounded 
%   below; set UB(i) = Inf if X(i) is unbounded above.
%
%   X = LSQLIN(C,d,A,b,Aeq,beq,LB,UB,X0) sets the starting point to X0.
%
%   X = LSQLIN(C,d,A,b,Aeq,beq,LB,UB,X0,OPTIONS) minimizes with the default
%   optimization parameters replaced by values in OPTIONS, an argument
%   created with the OPTIMOPTIONS function. See OPTIMOPTIONS for details.
%
%   X = LSQLIN(PROBLEM) solves the least squares problem defined in 
%   PROBLEM. PROBLEM is a structure with the matrix 'C' in PROBLEM.C, the 
%   vector 'd' in PROBLEM.d, the linear inequality constraints in 
%   PROBLEM.Aineq and PROBLEM.bineq, the linear equality constraints in 
%   PROBLEM.Aeq and PROBLEM.beq, the lower bounds in PROBLEM.lb, the upper 
%   bounds in PROBLEM.ub, the start point in PROBLEM.x0, the options 
%   structure in PROBLEM.options, and solver name 'lsqlin' in 
%   PROBLEM.solver. Use this syntax to solve at the command line a problem 
%   exported from OPTIMTOOL. 
%
%   [X,RESNORM] = LSQLIN(C,d,A,b) returns the value of the squared 2-norm 
%   of the residual: norm(C*X-d)^2.
%
%   [X,RESNORM,RESIDUAL] = LSQLIN(C,d,A,b) returns the residual: C*X-d.
%
%   [X,RESNORM,RESIDUAL,EXITFLAG] = LSQLIN(C,d,A,b) returns an EXITFLAG
%   that describes the exit condition. Possible values of EXITFLAG and the
%   corresponding exit conditions are
%  
%     1  LSQLIN converged to a solution X.
%     3  Change in the residual smaller that the specified tolerance.
%     0  Maximum number of iterations exceeded.
%    -2  Problem is infeasible.
%    -4  Ill-conditioning prevents further optimization.
%    -7  Magnitude of search direction became too small; no further 
%         progress can be made. The problem is ill-posed or badly 
%         conditioned.
%
%   [X,RESNORM,RESIDUAL,EXITFLAG,OUTPUT] = LSQLIN(C,d,A,b) returns a 
%   structure OUTPUT with the number of iterations taken in 
%   OUTPUT.iterations, the type of algorithm used in OUTPUT.algorithm, the 
%   number of conjugate gradient iterations (if used) in OUTPUT.cgiterations, 
%   a measure of first order optimality (large-scale algorithm only) in 
%   OUTPUT.firstorderopt, and the exit message in OUTPUT.message.
%
%   [X,RESNORM,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA] = LSQLIN(C,d,A,b) returns 
%   the set of Lagrangian multipliers LAMBDA, at the solution: 
%   LAMBDA.ineqlin for the linear inequalities C, LAMBDA.eqlin for the 
%   linear equalities Ceq, LAMBDA.lower for LB, and LAMBDA.upper for UB.
%
%   See also QUADPROG.

%   Copyright 1990-2015 The MathWorks, Inc.

defaultopt = struct( ...
    'Algorithm','trust-region-reflective', ...
    'Diagnostics','off', ...
    'Display','final', ...
    'JacobMult',[], ...
    'LargeScale','on', ...
    'MaxIter',200, ...
    'MaxPCGIter','max(1,floor(numberOfVariables/2))', ...
    'PrecondBandWidth',0, ...
    'TolCon',1e-8, ...
    'TolFun',100*eps, ...
    'TolFunValue',100*eps, ...
    'TolPCG',0.1, ...
    'TypicalX','ones(numberOfVariables,1)' ...
);


% If just 'defaults' passed in, return the default options in X
if nargin==1 && nargout <= 1 && isequal(C,'defaults')
   X = defaultopt;
   return
end

% Handle missing arguments
if nargin < 10
    options = [];
    if nargin < 9
        X0 = [];
        if nargin < 8
            ub = [];
            if nargin < 7
                lb = [];
                if nargin < 6
                    beq = [];
                    if nargin < 5
                        Aeq = [];
                        if nargin < 4
                            b = [];
                            if nargin < 3
                                A = [];
                            end
                        end
                    end
                end
            end
        end
    end
end

% Detect problem structure input
if nargin == 1
    if isa(C,'struct')
        [C,d,A,b,Aeq,beq,lb,ub,X0,options] = separateOptimStruct(C);
    else % Single input and non-structure.
        error(message('optimlib:lsqlin:InputArg'));
    end
end

% Prepare the options for the solver
[options, optionFeedback] = prepareOptionsForSolver(options, 'lsqlin');

if nargin == 0 
  error(message('optimlib:lsqlin:NotEnoughInputs'))
end
% Check for non-double inputs
msg = isoptimargdbl('LSQLIN', {'C','d','A','b','Aeq','beq','LB','UB','X0'}, ...
                                C,  d,  A,  b,  Aeq,  beq,  lb,  ub,  X0);
if ~isempty(msg)
    error('optimlib:lsqlin:NonDoubleInput',msg);
end


if nargout > 5
   computeLambda = true;
else 
   computeLambda = false;
end
if nargout > 4
   computeFirstOrderOpt = true;
else 
   computeFirstOrderOpt = false;
end

% Set up constant strings
activeSet =  'active-set';
trustRegion = 'trust-region-reflective'; 
unconstrained = 'mldivide';
interiorPoint = 'interior-point';
output.iterations = []; % initialize so that it will be the first field

% Options setup
largescale = strcmpi(optimget(options,'LargeScale',defaultopt,'fast'),'on');
% Read Algorithm
output.algorithm = optimget(options,'Algorithm',defaultopt,'fast');
if ~any(strcmpi(output.algorithm, {activeSet,trustRegion,interiorPoint}))
    error(message('optimlib:lsqlin:InvalidAlgorithm'));
end
algAndLargeScaleConflict = false;
if strcmpi(output.algorithm,trustRegion) && ~largescale
    % Conflicting options Algorithm='trust-region-reflective' and
    % LargeScale='off'. Choose active-set algorithm.
    % Warn later, not in case of early termination
    algAndLargeScaleConflict = true;
end
        
diagnostics = strcmpi(optimget(options,'Diagnostics',defaultopt,'fast'),'on');
mtxmpy = optimget(options,'JacobMult',defaultopt,'fast');
% Check if name clash
functionNameClashCheck('JacobMult',mtxmpy,'atamult','optimlib:lsqlin:JacobMultNameClash');

% Use internal Jacobian-multiply function if user does not provide JacobMult function 
if isempty(mtxmpy)
    mtxmpy = @atamult;
end

switch optimget(options,'Display',defaultopt,'fast')
case {'off','none'}
   verbosity = 0;
case {'iter','iter-detailed'}
   verbosity = 2;
case {'final','final-detailed'}
   verbosity = 1;
otherwise
   verbosity = 1;
end

% Set the constraints up: defaults and check size
[nineqcstr,numberOfVariables]=size(A);
neqcstr = size(Aeq,1);
ncstr = nineqcstr + neqcstr;

if isempty(C) || isempty(d)
   error(message('optimlib:lsqlin:FirstTwoArgsEmpty'))
else
   numberOfVariables = max([size(C,2),numberOfVariables]); % In case C is empty
end

[rows,cols]=size(C);
if length(d) ~= rows
   error(message('optimlib:lsqlin:InvalidCAndD'))
end

if length(b) ~= size(A,1)
   error(message('optimlib:lsqlin:InvalidAAndB'))
end

if length(beq) ~= size(Aeq,1)
   error(message('optimlib:lsqlin:InvalidAeqAndBeq'))
end

if ( ~isempty(A)) && (size(A,2) ~= cols) 
    error(message('optimlib:lsqlin:CAndA'))
end

if ( ~isempty(Aeq)) && (size(Aeq,2) ~= cols)
    error(message('optimlib:lsqlin:CAndAeq'))
end

if isempty(X0) && ~strcmpi(output.algorithm,interiorPoint)
    % This zero-valued X0 will potentially be changed in sllsbox or qpsub.
    % (This potentially temporary zero-valued x0 needed here for backwards 
    % compatibility because it's returned in output x if early termination 
    % occurs when bounds are inconsistent.)
    X0 = zeros(numberOfVariables,1);
    params.emptyInitialPoint = true;  % parameter passed to sllsbox    
else
    params.emptyInitialPoint = false; % parameter passed to sllsbox 
end
if isempty(A)
    A = zeros(0,numberOfVariables);
end
if isempty(b)
    b = zeros(0,1); 
end
if isempty(Aeq)
    Aeq = zeros(0,numberOfVariables);
end
if isempty(beq)
    beq = zeros(0,1);
end

% Set d, b and X to be column vectors
d = d(:);
b = b(:);
beq = beq(:);
X0 = X0(:);

[X0,lb,ub,msg] = checkbounds(X0,lb,ub,numberOfVariables);
if ~isempty(msg)
   exitflag = -2;
   [resnorm,residual,lambda]=deal([]);
   output.iterations = 0;
   output.algorithm = ''; % not known at this stage
   output.firstorderopt=[];
   output.cgiterations =[];
   output.message = msg;
   X=X0; 
   if verbosity > 0
      disp(msg)
   end
   return
end

hasLinearConstr = ~isempty(beq) || ~isempty(b);
% Test if C is all zeros or empty
if norm(C,'inf')==0 || isempty(C)
    if ~strcmpi(output.algorithm,interiorPoint)
        C = 0;
    else
        % C must be a sparse matrix of correct size to prevent an error in
        % the interior-point QP
        C = sparse(numel(d),numberOfVariables);
    end
end

caller = 'lsqlin';

% Test for constraints
if ~hasLinearConstr && all(isinf([lb;ub]))
    output.algorithm = unconstrained;
elseif strcmpi(output.algorithm,trustRegion)
    if algAndLargeScaleConflict
        warning(message('optimlib:lsqlin:AlgAndLargeScaleConflict'));
        output.algorithm = activeSet;
    elseif (rows < cols)
        warning(message('optimlib:lsqlin:MoreColsThanRows'));
        output.algorithm = activeSet;
    elseif hasLinearConstr
        warning(message('optimlib:lsqlin:LinConstraints'));
        output.algorithm = activeSet;
    end
end


if diagnostics
   % Do diagnostics on information so far
   gradflag = []; hessflag = []; constflag = false; gradconstflag = false; 
   non_eq=0;non_ineq=0; lin_eq=size(Aeq,1); lin_ineq=size(A,1); 
   XOUT=ones(numberOfVariables,1); funfcn{1} = []; confcn{1}=[];
   msg = diagnose('lsqlin',output,gradflag,hessflag,constflag,gradconstflag,...
      XOUT,non_eq,non_ineq,lin_eq,lin_ineq,lb,ub,funfcn,confcn);
end


switch output.algorithm;
case unconstrained
   % Disable the warnings about conditioning for singular and
   % nearly singular matrices
   warningstate1 = warning('off', 'MATLAB:nearlySingularMatrix');
   warningstate2 = warning('off', 'MATLAB:singularMatrix');
   X = C\d;
   % Restore the warning states to their original settings
   warning(warningstate1)
   warning(warningstate2)
   exitflag = 1;
   if computeFirstOrderOpt || computeLambda
       lambda.lower = [];
       lambda.upper = [];
       lambda.ineqlin = [];
       lambda.eqlin = [];
       output.iterations = 0;
       output.firstorderopt = [];
       output.constrviolation = [];
       output.message = '';
   end
case trustRegion
   params.verb = verbosity; % pack parameters together into struct
   defaultopt.TolFun = 100*eps;
   [X,resnorm,residual,firstorderopt,iterations,cgiterations,exitflag,lambda,msg]=...
      sllsbox(C,d,lb,ub,X0,params,options,defaultopt,mtxmpy,computeLambda,varargin{:});
   output.iterations = iterations;
   output.firstorderopt = firstorderopt;
   output.cgiterations = cgiterations;
   output.message = msg;
case activeSet   
   % Sparse inputs are not accepted by qpsub 
   if issparse(A) || issparse(C) || issparse(Aeq) 
      warning(message('optimlib:lsqlin:ConvertToFull'))
   end    
   % Create options structure for qpsub
   lsqlinoptions.MaxIter = optimget(options,'MaxIter',defaultopt,'fast');
   % A fixed constraint tolerance (eps) is used for constraint
   % satisfaction; no need to specify any value
   lsqlinoptions.TolCon = [];

   [X,lambdaqp,exitflag,output,~,~,msg]= ...
      qpsub(full(C),d,[full(Aeq);full(A)],[beq;b],lb,ub,X0,neqcstr, ...
            verbosity,caller,ncstr,numberOfVariables,lsqlinoptions);
   output.algorithm = activeSet; % qpsub overwrites output with output.iterations 
case interiorPoint
    defaultopt.MaxIter = 200;
    defaultopt.TolFun = 1e-8;
    defaultopt.TolCon = 1e-8;
    % Set ConvexCheck flag to notify QP solver that problem is convex.
    defaultopt.ConvexCheck = 0;
    % Need to compute product here since presolve requires the entire
    % matrix. It is also required to be sparse storage (even though it will
    % most likely be dense in terms of nonzeros) for internal computations.
    H = sparse(C'*C);
    f = -C'*d;
    A = sparse(A); Aeq = sparse(Aeq);
    % If the output structure is requested, we must reconstruct the
    % Lagrange multipliers in the postsolve. Therefore, set computeLambda
    % to true if the output structure is requested.
    flags.computeLambda = computeFirstOrderOpt; 
    flags.detailedExitMsg = false;
    flags.verbosity = verbosity;    

    [X,~,exitflag,output,lambda] = ipqpcommon(H,f,A,b,Aeq,beq,lb,ub,X0, ...
                                          flags,options,defaultopt,...
                                          optionFeedback,varargin{:});   
    output.algorithm = interiorPoint; % overwrite "interior-point-convex"
    output.cgiterations = [];
    if isempty(lambda)
        X = []; residual = []; resnorm = [];
        return;
    end
    
    % Presolve may have removed variables and constraints from the problem.
    % Postsolve will re-insert the primal and dual solutions after the main
    % algorithm has run. Therefore, constraint violation and first-order
    % optimality must be re-computed (below).                                  
end

if any(strcmpi(output.algorithm , {activeSet, interiorPoint, unconstrained}))
    residual = C*X-d;
    resnorm = sum(residual.*residual);
    
    if strcmpi(output.algorithm, activeSet)
        llb = length(lb);
        lub = length(ub);
        lambda.lower = zeros(llb,1);
        lambda.upper = zeros(lub,1);
        arglb = ~isinf(lb); lenarglb = nnz(arglb);
        argub = ~isinf(ub); lenargub = nnz(argub);
        lambda.eqlin = lambdaqp(1:neqcstr,1);
        lambda.ineqlin = lambdaqp(neqcstr+1:neqcstr+nineqcstr,1);
        lambda.lower(arglb) = lambdaqp(neqcstr+nineqcstr+1:neqcstr+nineqcstr+lenarglb);
        lambda.upper(argub) = lambdaqp(neqcstr+nineqcstr+lenarglb+1:neqcstr+nineqcstr+lenarglb+lenargub);
        output.message = msg;
        
        if exitflag == 1 % qpsub terminated successfully
            normalTerminationMsg = sprintf('Optimization terminated.');
            if verbosity > 0
                disp(normalTerminationMsg)
            end
            if isempty(msg)
                output.message = normalTerminationMsg;
            else
                % append normal termination msg to current output msg
                output.message = sprintf('%s\n%s',msg,normalTerminationMsg);
            end
        else
            output.message = msg;
        end
        % Form QP matrix if optimality is requested
        if computeFirstOrderOpt
            H = C'*C; f = -C'*d;
        end
    end
    % Compute constraint violation if the output structure is requested
    % Compute first order optimality if needed
    
    % For the interior-point algorithm, if no initial point was provided by
    % the user and the presolve has declared the problem infeasible or
    % unbounded, X will be empty. The lambda structure will also be empty,
    % so do not compute constraint violation or first-order optimality if
    % lambda is missing.
    if any(strcmpi(output.algorithm , {activeSet, interiorPoint})) && ...
            computeFirstOrderOpt
        output.constrviolation = norm([Aeq*X-beq;max([A*X - b;X - ub;lb - X],0)],Inf);        
        output.firstorderopt = computeKKTErrorForQPLP(H,f,A,b,Aeq,beq,lb,ub,lambda,X);
    end  
    output.cgiterations = [];
end
