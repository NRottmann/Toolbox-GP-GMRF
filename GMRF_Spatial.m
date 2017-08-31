function [mu,Sigma,time,opt_time,opt_p,opt_covParam] = GMRF_Spatial(p,x,s,y,varargin)

% Gaussian Process with build-in Gaussian Markov Random Field
%
% Syntax:
%   [mu,Sigma,time,opt_time,opt_p,opt_covParam] = GMRF_Spatial(p,x,s,y)
%   [mu,Sigma,time,opt_time,opt_p,opt_covParam] = GMRF_Spatial(p,x,s,y,'propertyname','propertyvalue',...)
%
% Description:
%   Generates a field estimation using a built-in Gaussian Markov Random
%   Field given by the generating points (p). To estimate the target values
%   at the wanted points (s) measurement values (y) at given points (x) are
%   needed. The ourput of the function is an estimation of the target
%   values at the points s. The mean and the variance for this points are
%   given as outputs.
%
% Input:
%   p: generating points, as a matrix with D x m (if generating points have
%   to be optimized, this will reduce to (2+2*D) x 1 vector with the input:
%   [D; m; min range dim 1; max range dim 1; min range dim 2; max range ... ])
%   x: measurement points as a matrix with D x n
%   s: unknown points as a matrix with D x l
%   y: measurements at points x with n x 1
% 
%   D: Dimension of the input space
%
% Propertyname/-value pairs:
%   noise - noise of the measurements (default: noise = 0.1)
%   CovFunc - name of the covariance function (string) which should be used
%   (default: se_kernel)
%   CovParam - Array of the Covariance Parameters, for further
%   information see Covariance Function description
%   PreFunc - name of the function to determine the precision matrix
%   (string) (default: radius_prec)
%   PreParam - array of the precision function parameters, for further
%   information see Precision Function description
%   improve - false for no improvement of parameter or generating points
%             true for improvement function (default: false)
%   optimthings: 
%       0: All parameter will be optimized (default)
%       1: Only the parameter of the covariance function will be optimized
%       2: Only the positions of the generating points will be optimized
%   method:
%       0: minimum search using fminunc with random start points
%       1: bayesian optimization with random start points
%       2: random search
%       3: genetic algortihms
%   maxTime: time which is used for the optimization
%   startParam: starting parameter for optimization process
%       default: 0 (random choose of starting parameters)
%       vector with size numHypParam x 1 (for ordering see code)
%   optimNumGen: default: false, if true the programm tries to find an
%   optimal number of generating points (only useable together with
%   generating points improvement, optimthings = 2) -> p has to be in the form
%           [min number; max number; min range x; max range x; min range y; max range y])
%   mean: the mean function which is used for prediction
%       zero: zero mean (default)
%       constant: a constant value as mean function produced due to the
%       mean value of the measurements
%   ComTime: choose which computation time should be given out
%           1: Computational time with generating of the covariance and
%           precision matrices (default)
%           2: Computational time only for the execution of the estimation
%           algorithm
%   
% Output:
%   mu - mean values for points s
%   Sigma - Variance and Covariance regarding the points s
%   time - needed computational time
%   opt_time - time needed for optimization
%   opt_p - optimized generating points
%   opt_covParam - optimized covariance parameters
%
% used subfunction: setargs
%
% Date: 24. August, 2016
% Author: Nils Rottmann

% Adding pathes
A = cell(4,1);
Full = mfilename('fullpath');
Name = mfilename;
base = Full(1:end-length(Name));
A{1} = [base 'subfunctions'];
A{2} = [base 'precision_functions'];
A{3} = [base 'covariance_functions'];
A{4} = [base 'optim_functions'];
for i=1:1:length(A)
    addpath(A{i});
end

% Default values
defaultargs = {'noise', 0.1, 'CovFunc', 'se_kernel', 'CovParam', [], ...
    'PreFunc', 'radius_prec', 'PreParam', [],'improve',false,...
    'optimthings',0,'maxTime',0.1,'method',0,'startParam',0,...
    'optimNumGen',false,'mean','zero','ComTime',1}; 
params = setargs(defaultargs, varargin);

% error checking
if params.optimthings ~= 1 && params.improve
    if size(x,1) ~= size(s,1)
        error('x and s must have the same number of rows')
    end
else
    if size(x,1) ~= size(s,1) || size(p,1) ~= size(s,1);
        error('p,x and s must have the same number of rows')
    end
end

if size(x,2) ~= size(y,1);
    error('The number of columns of x has to be the same as the numbe of rows in y')
end

% check for nans and delete them 
II = any(isnan(x),1);
III = isnan(y);
x(:,II) = [];
y(II)   = [];
x(:,III) = [];
y(III)   = [];

% Defining the call for the precision function and covariance function
Pre = str2func(params.PreFunc);
Cov = str2func(params.CovFunc);

%%%%%%%%%%%%%%%%%%%
% Optimizing parameters and generating points if wanted
%%%%%%%%%%%%%%%%%%%
% Getting the numbers of needed hyperparameter in the precision function
% and covariance function
test_vec_param = ones(length(x(:,1)),1);   % test vector for structure
if params.improve == true
    [~,hypCov] = Cov(test_vec_param,test_vec_param,'struct',true);
    if length(hypCov(1,:)) ~= 1 | iscell(hypCov)
        error('Wrong structure of default hyperparameters in used covariance function')
    end
    [~,hypPre] = Pre(test_vec_param,'struct',true);
    if length(hypPre(1,:)) ~= 1 | iscell(hypPre)
        error('Wrong structure of default hyperparameters in used precision function')
    end

    if isempty(params.PreParam)
        params.PreParam = hypPre;          % giving the default value
    end
    if isempty(params.CovParam)
        params.CovParam = hypCov;
    end
    
    % Checking if number of generating points should be optimized
    if params.optimthings == 2 && params.optimNumGen
        NumGen = p(1:2);
        FieldRange = p(3:6);
        [p_num_opt] = NumGenOpt(NumGen,x,y,params.noise,params.CovFunc,...
            params.CovParam,params.PreFunc,params.PreParam,FieldRange);
        disp('optimized number of generating points:')
        num_gen = length(p_num_opt(1,:))
        p = [length(p_num_opt(:,1)); num_gen; p(3:6)];
        params.startParam = zeros(2*num_gen,1);
        params.startParam(1:num_gen) = p_num_opt(1,:);
        params.startParam(num_gen+1:end) = p_num_opt(2,:);
    elseif params.optimNumGen
        error('No valid usage of optimNumGen. Optimthings = 2 and improve true is needed!')
    end
    
    
    % Optimizing
    [params.CovParam,p,opt_time] = hypParam_opt(p,x,y,params.noise,...
        params.CovFunc,params.CovParam,params.PreFunc,params.PreParam,...
        params.optimthings,params.method,params.maxTime,params.startParam);
    disp('Improved Covariance Function Parameters:')
    opt_covParam = params.CovParam
    disp('Improved positions generating points :')
    opt_p = p
    disp('Time needed for improvement :')
    opt_time
end

% rewrite parameters
sigma_n = params.noise;

%%%%%%%%%%%%%%%
% Getting the estimation
%%%%%%%%%%%%%%%
% start time measuring
if params.ComTime ~= 2
    tic
end

% Getting the mean function
if strcmp(params.mean,'zero')
    mu_func = @(z) (zeros(length(z(1,:)),1));
elseif strcmp(params.mean,'constant')
    mu_func = @(z) (mean(y) * ones(length(z(1,:)),1));
else
    error('Wrong input for variable mean!')
end

% Calculating the precision matrix
if isempty(params.PreParam)   % if no parameters are given take default values
    [Q] = Pre(p);
else
    [Q] = Pre(p,'PreParam',params.PreParam);
end

% creating covariance matrix between generating points and given points /
% unknown points
if isempty(params.CovParam)  % if no parameters are given take default values
    [G_x] = Cov(x,p);
    [G_s] = Cov(p,s);
else
    [G_x] = Cov(x,p,'CovParam',params.CovParam);
    [G_s] = Cov(p,s,'CovParam',params.CovParam);
end

% Algorithm to determine the mean value of the unknown points s
% start time measuring
if ComTime == 2
    tic
end

if sigma_n > 0.001
    Q_dach = Q + (sigma_n)^(-2) * (G_x') * G_x;
    y_dach = sigma_n^(-2) * G_x' * (y - mu_func(x));
else
    Q_dach = (G_x') * G_x;
    y_dach = G_x' * (y - mu_func(x));
    % Ensure positive definiteness
    Q_dach = 0.1*eye(length(Q_dach)) + Q_dach;  
end

L = chol(Q_dach,'lower');           % Cholesky Tranformation, L satisfies L*L' = Q_dach_inv
alpha = (L')\(L\y_dach);            % avoiding inverse transformation for computational effectiveness
mu = mu_func(s) + (G_s')*alpha;
beta = L \ (G_s);
Sigma = (beta')*beta;

time = toc;

end
