function [f] = get_f_bil(X,field,range)
% Getting the values at positions X (D x m) from field
% Field has to be given in a grid structure
% Bilinear Interpolation is used to determine the values https://en.wikipedia.org/wiki/Bilinear_interpolation
% Input:
% X: Positions (D x m)
% field: The field as a matrix (n_y x n_x)
% range: Range of the field as a vector (D x 1), for example if your field
% goes from 10m to 100m in x-direction your range would be 90m in
% x-direction
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
        x_real = 1 + X(1,i) * ((n_x-1)/range(1));           % exact positions in the matrix
        y_real = 1 + X(2,i) * ((n_y-1)/range(2));
        
        x = [floor(x_real); ceil(x_real)];                  % round (in both directions)
        y = [floor(y_real); ceil(y_real)];
        
        A = 1/((x(2) - x(1))*(y(2)- y(1)));            % interim values
        x_vec = [(x(2) - x_real) (x_real - x(1))];
        y_vec = [(y(2) - y_real); (y_real - y(1))];
        F = zeros(2,2);
        for ii=1:1:2
            for jj=1:1:2
                F(ii,jj) = field(y(jj),x(ii));
            end
        end            
        f(i) = A*x_vec*F*y_vec;                         % Bilinear Interpolation
    end
else
    error('Other dimensions not yet implemented!')
end


end

