function [mu,Sigma,Q_dach,y_dach] = GMRF_Sequential(Q_dach,y_dach,p,x,s,y,varargin)
% Gaussian Process with build-in Gaussian Markov Random Field implemented
% with zero mean and as a sequential algorithm
%
% Syntax:
%   [mu,Sigma,Q_dach,y_dach] = GMRF_Spatial(Q_dach,y_dach,p,x,s,y)
%   [mu,Sigma,Q_dach,y_dach] = GMRF_Spatial(Q_dach,y_dach,p,x,s,y,'propertyname','propertyvalue',...)
%
% Description:
%   Generates a field estimation using a built-in Gaussian Markov Random
%   Field given by the generating points (p). To estimate the target values
%   at the wanted points (s) measurement values (y) at given points (x) are
%   needed. The ourput of the function is an estimation of the target
%   values at the points s. The mean and the variance for this points are
%   given as outputs. Also Q_dach and y_dach to calculate the next time
%   step are given.
%
% Input:
%   Q_dach: in between precision matrix (if first no Q_dach available, set
%   Q_dach = [])
%   y_dach: in between output values (if first no y_dach available, set
%   y_dach = [])
%   p: generating points, as a matrix with D x m
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
%   
% Output:
%   mu - mean values for points s
%   Sigma - Variance and Covariance regarding the points s
%   Q_dach - in between precision function for next iteration
%   y_dach - in between output for next iteration
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
    'PreFunc', 'radius_prec', 'PreParam', []}; 
params = setargs(defaultargs, varargin);

% error checking
if size(x,1) ~= size(s,1) || size(p,1) ~= size(s,1);
    error('p,x and s must have the same number of rows')
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

% rewrite parameters
sigma_n = params.noise;

% checking if Q_dach and y_dach available and have the right size
num_p = length(p(1,:));
if isempty(Q_dach) && isempty(y_dach)
    disp('Calculate starting values for Q_dach and y_dach!')
    % Calculating the precision matrix
    if isempty(params.PreParam)   % if no parameters are given take default values
        [Q_dach] = Pre(p);
    else
        [Q_dach] = Pre(p,'PreParam',params.PreParam);
    end
    % initializing y_dach
    y_dach = zeros(num_p,1);
elseif (size(Q_dach) == [num_p num_p]) 
    if (size(y_dach) == [num_p 1])
        disp('Accpted Q_dach and y_dach!')
    else
        error('Wrong input of Q_dach and y_dach');
    end   
else
    error('Wrong input of Q_dach and y_dach');
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

% Updating precision matrix Q and vector y
Q_dach = Q_dach + (sigma_n)^(-2) * (G_x') * G_x;
y_dach = y_dach + sigma_n^(-2) * (G_x') * y;

L = chol(Q_dach,'lower');           % Cholesky Tranformation, L satisfies L*L' = Q_dach_inv
alpha = (L')\(L\y_dach);            % avoiding inverse transformation for computational effectiveness
mu = (G_s')*alpha;
beta = L \ (G_s);
Sigma = (beta')*beta;

end

