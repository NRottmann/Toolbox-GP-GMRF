function [K,varargout] = se_kernel(X1,X2,varargin)
% Squared-exponential covariance function (Kernel) with equal characteristic length
% in every dimension, it is a measure of similarity
%
% Formular:
%   k(x,x') = sigma_f^2 * exp(-(1/(2*(l^2))) * (x - x')^2) + sigma_n^2 kroneckerDelta_qp
%
% Syntax:
%   [K] = se_kernel(X1,X2)
%   [K,optionalOutputs] = se_kernel(X1,X2,'propertyname','propertyvalue',...)
%
% Description:
%   Generates a Covariance Matrix between the points given in X1 and X2.
%   This covariance matrix can be used for instance for Gaussian Processes
%
% Input: 
%   X1, X2: Matrices of data points with dimension D x n, D x m
%
% Propertyname/-value pairs:
%   CovParam - array of the covariance function parameters
%       2 x 1 array needed:
%       (1): scale factor, sigma_n (default: 1)
%       (2): characteristic length, l (default: 1)
%   struct - gives out the hyperparameters, additional output
%   argument needed (true or false, default: false)
%
% Output:
%   K - covariance matrix (n x m)
%
% Date: 17. August, 2016
% Author: Nils Rottman

% Default values
defaultargs = {'CovParam', [1; 1],'struct',0}; 
params = setargs(defaultargs, varargin);

% look if hyperparameters are wanted
if params.struct == true
    varargout{1} = params.CovParam;
end 

% error checking
if iscell(params.CovParam)
    error('covariance parameters are in a cell array')
end
if size(params.CovParam) ~= [2 1] 
    error('covariance parameters have not the right size')
end

% Rewriting parameters in variables
sigma_n = params.CovParam(1);
l = params.CovParam(2);

% computing the dimension of the kernel matrix
l1 = length(X1(1,:));
l2 = length(X2(1,:));

% providing the kernel matrix for effcient computing
K = zeros(l1,l2);

% Computing the kernel for each element of the kernel matrix
for i=1:1:l1
    for j=1:1:l2
        K(i,j) = sigma_n^2 * exp(-(1/(2*(l^2))) * (norm(X1(:,i)-X2(:,j)))^2);
    end
end
    
end

