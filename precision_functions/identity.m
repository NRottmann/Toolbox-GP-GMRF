function [Q,varargout] = identity(p,varargin)
% Precision function a Gaussian Process with build-in Gaussian Markov
% Random Field
%
% Formular
%           |N(i)| + c0     if j = i
% (Q)_ij =  -1              if j element of N(i)
%           0               otherwise
% distance function to define N(i): eps = {(i,j) | ||p_i - p_j|| < r}
%
% Syntax:
%       [Q] = radius_prec(p)
%       [Q,optionalOutputs] = radius_prec(p,'propertyname','propertyvalue',...)
%
% Description:
%   Generates the Precision Matrix for the Gaussian Process with build-in 
%   Gaussian Markov Random Field
%
% Input:
%   p: generating points , as a Matrix with D x m
%
% Propertyname/-value pairs:
%   PreParam - array of the precision function parameters
%       1 x 1 array: r -> characteristic radius
%   struct - gives out the hyperparameters, additional output
%   argument needed (true or false, default: false)
%
% Output:
%   Q - precision matrix (m x m)
%
% Date: 17. August, 2016
% Author: Nils Rottman

% Default values
defaultargs = {'PreParam', 0,'struct',0}; 
params = setargs(defaultargs, varargin);

% error checking
if iscell(params.PreParam)
    error('precision parameters are in a cell array')
end
if size(params.PreParam) ~= [1 1]
    error('Precision parameters have not the right size')
end

% look if hyperparameters are wanted
if params.struct == true
    varargout{1} = params.PreParam;
end 

l = length(p(1,:));
Q = eye(l);

end

