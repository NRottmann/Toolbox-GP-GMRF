function [Q,varargout] = pre_example(p,varargin)
% explain your precsion function
%
% equation
%    (Q)_ij = ...
%
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
%       (1): ...
%       (2): ...
%   struct - gives out the hyperparameters, additional output
%   argument needed (true or false, default: false)
%
% Output:
%   Q - precision matrix (m x m)
%
% Date: ...
% Author: ...

% Default values
defaultargs = {'PreParam', [...],'struct',0};
% add your default values
params = setargs(defaultargs, varargin);

% error checking
if iscell(params.PreParam)
    error('precision parameters are in a cell array')
end
if size(params.PreParam) ~= [... 1]
% add the number of cov.-param.
    error('Precision parameters have not the right size')
end

% look if hyperparameters are wanted
if params.struct == true
    varargout{1} = params.PreParam;
end 

% Rewriting parameters in variables
...
    
% Calculate Q
...



end

