function [CovParam_opt,p_opt,opt_time] = hypParam_opt(p,x,y,sigma_n,CovFunc,CovParam,PreFunc,PreParam,optimthings,method,maxTime,startParam)

% Hyperparameter and positioning of generating points optimization
% using the log likelihood method for the Gaussian Process with build-in 
% Gaussian Markov Random Field
%
% Syntax:
%   [CovParam_opt,p_opt] = hypParam_opt(p,x,y,sigma_n,CovFunc,PreFunc,optimthings)
%   [CovParam_opt,p_opt] = hypParam_opt(p,x,y,sigma_n,CovFunc,PreFunc,optimthings),'propertyname','propertyvalue',...)
%
% Description:
%   This function optimizes the hyperparameter for a Gaussian Process with
%   build-in GMRF. Therefore it uses the log likelihood which has to be
%   maximized. It is possible either to optimize just the covariance
%   function parameters, just the positions of the generating points or
%   both.
%
% Input:
%   p: generating points, as a matrix with D x m (if generating points have
%   to be optimized, this will reduce to (2+2*D) x 1 vector with the input:
%   [D; m; min range dim 1; max range dim 1; min range dim 2; max range ... ])
%   x: measurement points as a matrix with D x n
%   y: measurements at points x with n x 1
%   sigma_n: noise factor in the data
%   CovFunc: name of the covariance function to determine the covariance matrix (see
%            folder)
%   CovParam: Hyperparameter for the covariance function
%   PreFunc: name of the precision function to determine the precision
%            matrix
%   PreParam: Hyperparameter for the Precision function
%   optimthings: 
%       0: All parameter will be optimized
%       1: Only the parameter of the covariance function will be optimized
%       2: Only the positions of the generating pointd will be optimized
%   maxTime: maximum time for optimization
%   startParam: defines the start parameter for the estimation process
%   
% Output:
%   CovParam_opt: optimized Hyperparameter
%   p_opt: optimized positions of the generating points

% rewriting p iff generating points positions have to be optimized
if optimthings == 0 || optimthings == 2
    D = p(1);
    vec = size(p);
    if vec(1) ~= (2+2*D) || vec(2) ~= 1
        error('wrong size of the precision vector for generating points optimization!')
    end
    m = p(2);
    num_gen = D*m;
    for f=1:1:D
        range(:,f) = p((3+(f-1)*2):(4+(f-1)*2));
    end
else
    range = 10;
end

% Calculating the number of hyperparameters and defining the objective
% function
num_cov = length(CovParam);
if optimthings == 0
    num_all = num_gen + num_cov;
    objFunc = @(obj)logLikelihood_GMRF(obj(num_cov+1:num_all),x,y,sigma_n,CovFunc,obj(1:num_cov),PreFunc,PreParam,D,m,optimthings);
    objFunc_bay = @(obj)logLikelihood_GMRF_bay(obj(num_cov+1:num_all),x,y,sigma_n,CovFunc,obj(1:num_cov),PreFunc,PreParam,D,m);
    objFunc_gen = @(obj)logLikelihood_GMRF(obj(num_cov+1:num_all),x,y,sigma_n,CovFunc,obj(1:num_cov)',PreFunc,PreParam,D,m,optimthings);
elseif optimthings == 1
    num_all = num_cov;
    [D, m] = size(p);
    for c=1:1:D
        p_ord((1+(m*(c-1))):(c*m)) = p(c,:);
    end
    objFunc = @(obj)logLikelihood_GMRF(p_ord,x,y,sigma_n,CovFunc,obj,PreFunc,PreParam,D,m,optimthings);
    objFunc_bay = @(obj)logLikelihood_GMRF_bay(p_ord,x,y,sigma_n,CovFunc,obj,PreFunc,PreParam,D,m);
    objFunc_gen = @(obj)logLikelihood_GMRF(p_ord,x,y,sigma_n,CovFunc,obj',PreFunc,PreParam,D,m,optimthings);
elseif optimthings == 2
    num_all = num_gen;
    objFunc = @(obj)logLikelihood_GMRF(obj,x,y,sigma_n,CovFunc,CovParam,PreFunc,PreParam,D,m,optimthings);
    objFunc_bay = @(obj)logLikelihood_GMRF_bay(obj,x,y,sigma_n,CovFunc,CovParam,PreFunc,PreParam,D,m);
    objFunc_gen = @(obj)logLikelihood_GMRF(obj,x,y,sigma_n,CovFunc,CovParam,PreFunc,PreParam,D,m,optimthings);
end

%%%%%%%%%%%%%%%%%
% Doing the optimization
%%%%%%%%%%%%%%%%%
if method == 0          % random points with minimum search using fminunc
    disp('start optimization using fminunc')
    % options
    opts.Algorithm = 'trust-region';
    opts.GradObj = 'on';
    % decision what has to be optimized
    opt_time = 0;
    fval_old = inf;
    obj_new = zeros(num_all,1);
    while opt_time < maxTime
        tic;
        if startParam == 0
            [obj_0] = random(num_all,num_cov,D,m,range,optimthings);
        else
            obj_0 = startParam;
        end
        [obj_new,fval_new] = fminunc(objFunc,obj_0,opts);
        if fval_new < fval_old
            fval_old = fval_new;
            obj_old = obj_new;
        end
        dtime = toc;
        opt_time = opt_time + dtime
    end
    disp('end optimization using fminunc')
      
elseif method == 1      % bayesian optimization
    disp('start optimization using bayesian optimization')
    [x_bay_new] = random(num_all,num_cov,D,m,range,optimthings);
    x_bay = x_bay_new;
    y_bay = objFunc_bay(x_bay);

    % Starting the iteration loop
    opt_time = 0;
    while opt_time < maxTime
        tic;
        
        if startParam == 0
            [x_bay_rand] = random(num_all,num_cov,D,m,range,optimthings);
        else
            error('no start parameter allowed using this optimization method')
        end
        
        objFunc_GP = @(obj)BayOpt_objFun(x_bay,obj,y_bay,num_all);      % Defining the objective function
        opts_bay = optimoptions(@fminunc,'Algorithm','quasi-newton','GradObj','on');
        [x_bay_new,~] = fminunc(objFunc_GP,x_bay_rand,opts_bay);
        y_bay_new = objFunc_bay(x_bay_new);
        x_bay = [x_bay x_bay_new];
        y_bay = [y_bay; y_bay_new];

        dtime = toc;
        opt_time = opt_time + dtime
    end
    
    [~,ind] = min(y_bay);           % find minimum value in y
    obj_old = x_bay(:,ind);
    disp('end optimization using bayesian optimization')
    
elseif method == 2          % random tryouts
    disp('start optimization using random tryouts')
    opt_time = 0;
    fval_old = inf;
    while opt_time < maxTime
        tic;
        if startParam == 0
            [obj_new] = random(num_all,num_cov,D,m,range,optimthings);
        else
            error('no start parameter allowed using this optimization method')
        end
        [fval_new] = objFunc(obj_new);
        if fval_new < fval_old
            fval_old = fval_new;
            obj_old = obj_new;
        end
        dtime = toc;
        opt_time = opt_time + dtime
    end
    disp('end optimization using random tryouts')
    
elseif method == 3      % Genetic algorithms
    % Creating upper and lower bound for search
    LB = zeros(num_all,1);
    UB = zeros(num_all,1);
    if optimthings == 0        
        LB(1:num_cov) = 0;      % defining the upper and lower bound for covariance hyperparameters
        UB(1:num_cov) = 10;
        for ll = 1:1:D
            for jj=1:1:m
                LB(num_cov+(ll-1)*m+jj) = range(1,ll);
                UB(num_cov+(ll-1)*m+jj) = range(2,ll);
            end
        end
    elseif optimthings == 1
        LB(1:num_cov) = 0;      % defining the upper and lower bound for covariance hyperparameters
        UB(1:num_cov) = 10;
    elseif optimthings == 2
        for ll = 1:1:D
            for jj=1:1:m
                LB((ll-1)*m+jj) = range(1,ll);
                UB((ll-1)*m+jj) = range(2,ll);
            end
        end
    else
        error('Wrong input for optimthings!')
    end
    
    nvars = double(num_all);
    ga_options = gaoptimset('TimeLimit',maxTime);
    obj_gen = ga(objFunc_gen,nvars,[],[],[],[],LB,UB,[],ga_options);
    obj_old = obj_gen';
    
    opt_time = maxTime;
    
else
    error('no valid value for method!')
end




% Redistributing the optimized parameters
if optimthings == 0
    CovParam_opt = obj_old(1:num_cov);
    for l=1:1:D
        p_opt(l,:) = obj_old((num_cov+1+(m*(l-1))):(num_cov+l*m));
    end
elseif optimthings == 1
    CovParam_opt = obj_old;
    p_opt = p;
elseif optimthings == 2
    for l=1:1:D
        p_opt(l,:) = obj_old((1+(m*(l-1))):(l*m));
    end
    CovParam_opt = CovParam;
end



end
    



