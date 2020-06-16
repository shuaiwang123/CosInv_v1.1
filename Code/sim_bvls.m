%
% x=bvls(A,b,l,u)
%
% Solves a bounded variables least squares problem
%
%   min || Ax-b ||
%          l <= x <= u
%
% The algorithm is based on the bounded variables least squares algorithm
% of Stark and Parker.  See
%
%  P.B. Stark and R. L. Parker, "Bounded-Variable Least-Squares: An Algorithm
%  and Applications", Computational Statistics 10:129-141, 1995.
%%
% introduced into PSOKINV by Feng,W.P., @ Glasgow, 2012-09-09
%
%%
%
function x = sim_bvls(A,b,l,u,opt)
%
% Get the size of A for future reference.
%
%[m,n]=size(A,2);
n = size(A,2);
%
% Initialize a bunch of variables.  This speeds things up by avoiding
% realloation of storage each time an entry is added to the end of a vector.
%
oopslist = [];
state    = zeros(n,1);
x        = zeros(n,1);
atbound  = [];
between  = [];
criti    = 0;
myeps    = 1.0e-15;
%
dims     = size(A);
maxiter  = 10*dims(2);       % add by W.P. Feng, 2012-09-09
                       % change back to 1000 from 50000, it will take long time if the input is a
                       % large data set... @ GU, 2013-03-24
                       %
%
if maxiter > 10000
    maxiter = 1000;
end
%
%
if nargin > 4
    if(isfield(opt,'maxiter'))
        maxiter    = opt.maxiter;
    end
    if(isfield(opt,'myeps'))
        myeps = opt.myeps;
    end
    if(isfield(opt,'x0'))
        x = opt.x0;
    end
end
%
% setup an initial solution with all vars locked at a lower or upper bound.
% set any free vars to 0.
%
for i=1:n
    if ((u(i) >= +Inf) && (l(i) <= -Inf))
        x(i)     = 0;
        state(i) = 0;
        between  = [between i];
        continue;
    end
    if (u(i) >= +Inf)
        x(i)=l(i);
        state(i)=1;
        atbound=[atbound i];
        continue;
    end
    if (l(i) <= -Inf)
        x(i)=u(i);
        state(i)=2;
        atbound = [atbound i];
        continue;
    end
    if (abs(l(i)) <= abs(u(i)))
        x(i)=l(i);
        state(i)=1;
        atbound=[atbound i];
        continue;
    else
        x(i)=u(i);
        state(i)=2;
        atbound=[atbound i];
        continue;
    end
end
%
% Print out some info.
%
%  atbound
%  between
%  x
%  pause;

%
% The main loop.  Stop after 10*n iterations if nothing else.
%
iter = 0;
while (iter < maxiter)
    %
    iter = iter+1;
    %
    % optimality test.
    %
    grad = A'*(A*x-b);
    %
    % We ignore any variable that has already been tried and failed.
    %
    grad(oopslist) = zeros(size(oopslist));
    %
    % Check for optimality.
    %
    done = 1;
    for i=1:n
        if ((abs(grad(i)) > (1+norm(b))*myeps) && (state(i)==0))
            done=0;
            break;
        end
        if ((grad(i) < 0) && (state(i)==1))
            done=0;
            break;
        end
        if ((grad(i) > 0) && (state(i)==2))
            done=0;
            break;
        end
    end
    
    if (done == 1)
        return;
    end
    %
    % Not optimal, so we need to free up some vars at bounds.
    %
    newi = 0;
    newg = 0.0;
    for i=atbound
        %
        % Don't free up a variable that was just locked at its bound.
        %
        if (i==criti)
            continue;
        end
        %
        % Look for the locked variable with the biggest gradient.
        %
        if ((grad(i) > 0) && (state(i)==2))
            if (abs(grad(i))>newg)
                newi=i;
                newg=abs(grad(i));
            end
        end
        if ((grad(i) < 0) && (state(i)==1))
            if (abs(grad(i))>newg)
                newi=i;
                newg=abs(grad(i));
            end
        end
    end
    
    %
    % Free the locked variable with the biggest gradient if there is one.
    %
    if (newi ~= 0)
        atbound    = remove(atbound,newi);
        state(newi)= 0;
        between    = [between newi];
    end
    %
    % Make sure the projected problem is nontrivial.
    %
    if (isempty(between))
        disp('Empty projected problem');
        continue;
    end
    %
    % Construct the new projected problem.
    %
    Aproj = A(:,between);
    An    = A(:,atbound);
    if (~isempty(atbound))
        bproj = b-An*x(atbound);
    else
        bproj = b;
    end
    %
    % Solve the projected problem.
    %
    z      = Aproj\bproj;
    xnew   = zeros(size(x));
    xnew(atbound) = x(atbound);
    xnew(between) = z;
    
    if ((newi ~= 0) && (((xnew(newi)<= l(newi)) && (x(newi)==l(newi))) || ...
            ((xnew(newi) >= u(newi)) && (x(newi)==u(newi)))))
        %
        % Ooops- the freed variable wants to go beyond its bound.  Lock it down
        % and look for some other variable.
        %
        oopslist = [oopslist; newi];
        if ((xnew(newi) <= l(newi)) && (state(newi)==1))
            state(newi)=1;
            x(newi)=l(newi);
        end
        if ((xnew(newi) >= u(newi)) && (state(newi)==2))
            state(newi)=2;
            x(newi)=u(newi);
        end
        atbound = [atbound newi];
        between = remove(between,newi);
        continue;
    end
    %
    % We've got a good variable freed up.  Reset the oopslist.
    %
    oopslist=[];
    %
    % Move as far as possible towards the optimal solution to the projected
    % problem.
    %
    alpha = 1;
    for i = between
        if (xnew(i) > u(i))
            newalpha=min([alpha,(u(i)-x(i))/(xnew(i)-x(i))]);
            if (newalpha < alpha)
                criti=i;
                crits=2;
                alpha=newalpha;
            end
        end
        if (xnew(i) < l(i))
            newalpha=min([alpha,(l(i)-x(i))/(xnew(i)-x(i))]);
            if (newalpha < alpha)
                criti = i;
                crits = 1;
                alpha=newalpha;
            end
        end
    end
    
    %  alpha
    %
    % Take the step.
    %
    x = x + alpha*(xnew-x);
    %
    % Update the state of variables.
    %
    if (alpha < 1)
        between = remove(between,criti);
        atbound = [atbound criti];
        state(criti)=crits;
    end
    
    for i = 1:n
        if (x(i) >= u(i))
            x(i)     = u(i);
            state(i) = 2;
            if ((isempty(between)) || (~isempty(find(between==i, 1))))
                between = remove(between,i);
            end
            if (isempty(find(atbound==i, 1)))
                atbound = [atbound i];
            end
        end
        if (x(i) <= l(i))
            x(i)     = l(i);
            state(i) = 1;
            if ((isempty(between)) || (~isempty(find(between==i, 1))))
                between = remove(between,i);
            end
            if (isempty(find(atbound==i, 1)))
                atbound = [atbound i];
            end
        end
    end
    
    %
    %
    % Go back and do it again.
    %
end
%disp('BVLS Exceeded max iters')

