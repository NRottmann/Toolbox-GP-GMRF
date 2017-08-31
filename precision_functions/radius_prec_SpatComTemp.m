function [Q,varargout] = radius_prec_SpatComTemp(p,varargin)
% Precision function a Gaussian Process with build-in Gaussian Markov
% Random Field
%
% Formular
%           |N(i)| + c0     if j = i
% (Q)_ij =  -1              if j element of N(i)
%           0               otherwise
% distance function to define N(i): eps = {(i,j) | ||p_i - p_j|| < r and the same time step or one previous}
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
defaultargs = {'PreParam', 3,'struct',0}; 
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

% Rewriting parameters in variables
r = params.PreParam;

% Factor for stabilizing positive definitness
c0 = 0.1; 
% number of generating points
l = length(p(1,:));     
% Calculation of the degree of node i for every generating point
N = false(l); % creates a lxl matrix with logical false (0) in it
% N: true or false if j is element of N(i)
N_i = zeros(l,1);
for i = 1:1:l
    deg = 0;
    for j = 1:1:l
        if norm(p(1:2,i) - p(1:2,j)) < r && (p(3,i) - p(3,j)) < 1.1
            N(i,j) = true;
            deg = deg + 1;
        end
    end
    N_i(i) = deg;
end

% Calculating the precision matrix
l_p = length(p(1,:));
Q = zeros(l_p,l_p);
for i = 1:1:l_p
    for j = 1:1:l_p
        if j == i
            Q(i,j) = N_i(i) + c0;
        elseif N(i,j) == true
            Q(i,j) = -1;
        else
            Q(i,j) = 0;
        end
    end
end

end

