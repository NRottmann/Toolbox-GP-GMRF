function [x,m] = MeshGen(n,range)
% Generating Mesh Points for the GMRF, maximum Dimension is D = 3
% Info
% D: Dimension of the field 
% Input:
% n: Number of mesh points in each direction (mesh) (D x 1)
% range: range of the field in each direction (D x 1)
% Output:
% x: mesh points in a D x m
% m: total number of mesh points

% Checking if Dimensions agree
D = length(n);
if D == length(range) && length(n(1,:)) == 1 && length(range(1,:)) == 1 && D < 4
    disp('Dimensions of input vectors checked!')
else
    error('Dimensions of input vectors are not correct!')
end

% mesh points in each axis with equal distance in between
X = cell(D,1);
for i=1:1:D
    X{i} = linspace(0,range(i),n(i));
end
disp('points on axis created!')

% Combining all input variables to one structure to create the mesh
% points
m = prod(n);    
x = zeros(D,m);
if D == 1
    x = X{1};
elseif D == 2
    for i=1:1:n(1)
        zw = ones(1,n(2))*X{1}(i);
        x(:,(1+(n(2)*(i-1))):(n(2)*i)) = [zw; X{2}];
    end
elseif D == 3
    error('Not implemented yet!')
end
    
end

