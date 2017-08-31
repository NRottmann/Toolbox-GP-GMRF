function [mu,Sigma] = GP_opt(x,s,y)
% Gaussian Processes for Bayesian Optimization
% This is a optimized version of Gaussian Processes for the usage in
% Bayesian Optimization. Call with hypParam_opt.m

% Covariance Matrices
[K] = se_kernel(x,x);
[K_s] = se_kernel(x,s);
[K_ss] = se_kernel(s,s);

% Ensure positive definiteness
K = 0.1*eye(length(K)) + K;                    

% Calculating mean vector and covariance matrix
I = eye(length(K));
L = chol(K,'lower');    % Cholesky Tranformation, L satisfies L*L' = K
alpha = (L')\(L\y);                     % avoiding inverse transformation for computational effectiveness
mu = (K_s')*alpha;
beta = L \ K_s;
Sigma = K_ss - (beta')*beta;


end

