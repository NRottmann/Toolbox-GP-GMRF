function [gamma] = variogram(M,k,h1,h2,fsize)
% Calculates the directional semi-variogram for a 2D field in which the
% Data is Aligned and Regularly Spaced. The field has to be given as a
% Matrix.
%
% Syntax:
%
%
% Input:
%   M: Matrix with field, size (n_y,n_x)
%   k: Defines the the maximum k for k*h_i by calculating the variogram
%   h1: distance in x-direction between two neighbouring points
%   h2: distance in y-direction between two neighbouring points

% Get the length of the field
lx = length(M(1,:));
ly = length(M(:,1));

% Calculate the h1, h2, h3, h4 for directions alpha1, alpha2, alpha3 and
% alpha4
h3 = sqrt(h1^2 + h2^2);
h4 = h3;

% calculate the alpha [degree]
alpha1 = 0;
alpha2 = -90;
alpha3 = 135;
alpha4 = 45;

% Initialize the variogram variables
gamma = cell(4,1);
for i=1:1:4
    gamma{i} = zeros(k+1,1);
end


% Compute the variogram values
% direction alpha1
for ad=1:1:k
    sum = 0;
    n = 0;
    for j=1:1:lx
        for i=1:1:ly
            if (j+ad) > lx || isnan(M(i,j+ad)) || isnan(M(i,j))
            else
                sum = sum + (M(i,j+ad) - M(i,j))^2;
                n = n + 1;
            end
        end
    end
    gamma{1}(ad+1) = 1/(2*n) * sum;
end

% direction alpha2
for ad=1:1:k
    sum = 0;
    n = 0;
    for j=1:1:lx
        for i=1:1:ly
            if (i+ad) > ly || isnan(M(i+ad,j)) || isnan(M(i,j))
            else
                sum = sum + (M(i+ad,j) - M(i,j))^2;
                n = n + 1;
            end
        end
    end
    gamma{2}(ad+1) = 1/(2*n) * sum;
end

% direction alpha3
for ad=1:1:k
    sum = 0;
    n = 0;
    for j=1:1:lx
        for i=1:1:ly
            if (i-ad) < 1 || (j-ad) < 1 || isnan(M(i-ad,j-ad)) || isnan(M(i,j))
            else
                sum = sum + (M(i-ad,j-ad) - M(i,j))^2;
                n = n + 1;
            end
        end
    end
    gamma{3}(ad+1) = 1/(2*n) * sum;
end

% direction alpha4
for ad=1:1:k
    sum = 0;
    n = 0;
    for j=1:1:lx
        for i=1:1:ly
            if (i-ad) < 1 || (j+ad) > lx || isnan(M(i-ad,j+ad)) || isnan(M(i,j))
            else
                sum = sum + (M(i-ad,j+ad) - M(i,j))^2;
                n = n + 1;
            end
        end
    end
    gamma{4}(ad+1) = 1/(2*n) * sum;
end

% generate x-vectors for plotting
x1 = 0:h1:k*h1;
x2 = 0:h2:k*h2;
x3 = 0:h3:k*h3;
x4 = x3;

% plot the results
% Normal Plot
figure
plot(x1,gamma{1},x2,gamma{2},x3,gamma{3},x4,gamma{4})
legend('\alpha_1','\alpha_2','\alpha_3','\alpha_4')
xlabel('|h|')
ylabel('\gamma')
title('Semi-Variogram')
set(gca,'FontSize',fsize)

% Contour Plot
X = [x1'; zeros(length(x2),1); -x1'; x1'; -x1'];
Y = [zeros(length(x1),1); x2'; x2'; x2'; zeros(length(x1),1)];
Z = [gamma{1}; gamma{2}; gamma{3}; gamma{4}; gamma{1}];
xnodes = -k*h1:h1:k*h1;
ynodes = 0:h2:k*h2;
[zg,xg,yg] = gridfit(X,Y,Z,xnodes,ynodes);

zg_add = rot90(zg,2);
zg_add(1,:) = [];
zg = [zg_add; zg];

yg_add = -flipud(yg(:,1));
yg_add(end) = [];
yg = [yg_add; yg(:,1)];

xg = xg(1,:);

max_value = max([max(gamma{1}),max(gamma{2}),max(gamma{3}),max(gamma{4})]);

figure
contourf(xg,yg,zg,30,'LineStyle','none')
caxis([0 max_value])
colormap(cool)
colorbar
title('Directional Semi-Variogram')
xlabel('h_x')
ylabel('h_y')
hold on
plot(x1',zeros(length(x1),1),'--k',zeros(length(x2),1),-x2','--k',...
    -x1',x2','--k',x1',x2','--k')
set(gca,'FontSize',fsize)




end

