function [f] = get_f_NN(X,field,range)
% Getting the values at positions X (D x m) from field
% Gets the nearest neighbour (NN) to the given points
% Input:
% X: Positions (D x m)
% field: The field as a matrix (n_y x n_x)
% range: Range of the field as a vector (D x 1)
% Output:
% f: values according to the points in X (D x m)

n = size(field);
n_x = n(2);
n_y = n(1);
m = length(X(1,:));
D = length(X(:,1));

if D == length(range) && length(range(1,:)) == 1
else
    error('Dimensions disagree!')
end

f = zeros(m,1);
if D == 2
    for i=1:1:m
        x_real = 1 + X(1,i) * ((n_x-1)/range(1));           % excat positions in the matrix
        y_real = 1 + X(2,i) * ((n_y-1)/range(2));
             
        f(i) = field(round(y_real),round(x_real));                         % Bilinear Interpolation
    end
else
    error('Other dimensions not yet implemented!')
end


end

