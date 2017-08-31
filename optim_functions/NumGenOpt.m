function [p] = NumGenOpt(NumGen,x,y,sigma_n,CovFunc,CovParam,PreFunc,PreParam,FieldRange)
% 
% Optimization of the number of Generating Points used in the GP with 
% built-in GMRF using an Information Criterion (here Akaike)
% AIC = 2k - 2ln(L)  
% L: Likelihood function; k: number of generating points
% Only possible for Generating points in 2D


% Input:
% NumGen: Number of Generating Points which should be tested, the values
% are guiding values for the function [NumMin; NumMax];
% FieldRange: Range of the field [min range x; max range x; min range y; max range y])
% x: given points as a matrix with D x n
% f: given values for points x
% sigma_n: noise factor in the data
% CovFunc: name of the function to determine the covariance matrix (see
% folder)
% HParam = Hyperparamter as a vector, see chosen covarianc function
% PrecFunc: name of the function to determine the precision matrix (see
% folder)
% PreParam = Input paramter as a cell array, see chosen precision function
% FieldRange: vector (D:1) which gives the range of the field

% Output:
% p: Generating Points, as a matrix with D x m


% test if the range of generating points is correctly given
if length(NumGen(:,1)) == 2 && length(NumGen(1,:)) == 1
    if NumGen(1) < NumGen(2)
        disp('Starting optimization of the number of gen. points!')
    else
        error('The minimum number and maximum number of generating points are given are in correct given!')
    end
else
    error('The minimum number and maximum number of generating points are given are in correct given!')
end
if size(FieldRange) == [4,1]
    if FieldRange(1) < FieldRange(2) && FieldRange(3) < FieldRange(4)
    else
        error('Range of the field is incorrect!')
    end     
else
    error('Range of the field is incorrect!')
end

% Define variables
l01 = FieldRange(1);
l02 = FieldRange(3);
l1 = FieldRange(2) - FieldRange(1);
l2 = FieldRange(4) - FieldRange(3);

% Starting the iteration to optimize the number of gen points
N = NumGen(1);                      % Number of wanted generating points
iter = NumGen(2) - N;
N_real = 0;
AIC = inf;
for i=1:1:iter
    n = round(sqrt(l1*N/l2));       % number of points in every line in x-direction
    m = round(sqrt(l2*N/l1));       % number of points in every line in y-direction
    
    N_real_new = n*m;
    
    if N_real_new > N_real
        % Creating the gen. points
        p_new = zeros(2,N_real_new);
        dist1 = l1/(n+1);
        dist2 = l2/(m+1);
        for j=1:1:n
            for k=1:1:m
                p_new(:,((j-1)*m + k)) = [l01 + j*dist1; l02 + k*dist2];   % coordinates
            end
        end

        % define functions to call them
        Pre = str2func(PreFunc);
        Cov = str2func(CovFunc);

        % Calculation of the COvariance Matrix and the precision matrix
        [Lambda] = Cov(x,p_new,'CovParam',CovParam);
        [Q] = Pre(p_new,'PreParam',PreParam);

        % Calculating the Covariance Matrix for y
        C = Lambda*Q^(-1)*Lambda' + sigma_n*eye(length(Lambda(:,1)));
        C_inv = C^(-1);

        % Calculating the marginal likelihood multiplied by -1
        n = length(y);
        logLik = - (0.5)*y'*C_inv*y - 0.5*log(norm(C)) - (n/2)*log(2*pi); 

        % Calculating the Akaike Value
        AIC_new = 2*m - 2*logLik;

        if AIC > AIC_new
            AIC = AIC_new;
            p = p_new;
        end
    end
    N_real = N_real_new;
    N = N + 1;
end

