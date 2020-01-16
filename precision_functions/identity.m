function [Q,varargout] = identity(p,varargin)
% Precision function a Gaussian Process with build-in Gaussian Markov
% Random Field, this is simply the identity matrix
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

