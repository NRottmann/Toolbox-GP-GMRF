function [K,varargout] = distr_approach(X1,X2,varargin)
% distributed covariance function      % explain your covariance function
%
% equation:
%   k(x,x') = (1 - h)cos(pi*h) + 1/pi sin(pi*h), h <= 1
%                  0                               otherwise
%   h = ||x - x'|| / r
%
% Syntax:
%   [K] = se_kernel(X1,X2)
%   [K,optionalOutputs] = se_kernel(X1,X2,...
%               'propertyname','propertyvalue',...)
%
% Description:
%   Generates a Covariance Matrix between the points given in X1 and X2.
%   This covariance matrix can be used for instance for Gaussian Processes.
%
% Input: 
%   X1, X2: Matrices of data points with dimension D x n, D x m
%
% Propertyname/-value pairs:
%   CovParam - array of the covariance function parameters
%       1 x 1 array needed:
%       (1): r
%   struct - gives out the hyperparameters, additional output
%   argument needed (true or false, default: false) (this is needed for
%   optimization scheme
%
% Output:
%   K - covariance matrix (n x m)
%
% Date: 01.11.16
% Author: Nils Rottmann

% Default values
defaultargs = {'CovParam', [1],'struct',0};   % add your default values
params = setargs(defaultargs, varargin);

% look if hyperparameters are wanted
if params.struct == true
    varargout{1} = params.CovParam;
end 

% error checking
if iscell(params.CovParam)
    error('covariance parameters are in a cell array')
end
if isscalar(params.CovParam)     % add the number of hyp.-param.
else
    error('covariance parameters have not the right size')
end

% Rewriting parameters in variables
r =  params.CovParam;  % here you can rename your hyperparameter if wanted

% computing the dimension of the kernel matrix
l1 = length(X1(1,:));
l2 = length(X2(1,:));

% providing the kernel matrix for effcient computing
K = zeros(l1,l2);

% Computing the kernel for each element of the kernel matrix
for i=1:1:l1
    for j=1:1:l2
        h = norm(X1(:,i) - X2(:,j))/r;
        if h > 1
            K(i,j) = 0;    
        else
            K(i,j) = (1 - h)*cos(pi*h) + (1/pi) * sin(pi*h);
        end
    end
end
    
end

