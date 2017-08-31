function [f] = logLikelihood_GMRF_bay(p_ord,x,y,sigma_n,CovFunc,CovParam,PreFunc,PreParam,D,m)
    % Loglikelihood Function for GP with built-in GMRF using Bayesian
    % Optimization
    % This is the Loglikelihood function which has to be minimized for
    % finding optimal parameter for a GP with built-in GMRF if Bayesian 
    % Optimization is used. Call due to hypParam_opt.m
    
    % define functions to call them
    Pre = str2func(PreFunc);
    Cov = str2func(CovFunc);

    % generating points have a predefined ordering
    for c=1:1:D
        p(c,:) = p_ord((1+(m*(c-1))):(c*m));
    end
    
    % Calculation of the COvariance Matrix and the precision matrix
    [Lambda] = Cov(x,p,'CovParam',CovParam);
    [Q] = Pre(p,'PreParam',PreParam);

    % Calculating the Covariance Matrix for y
    C = Lambda*Q^(-1)*Lambda' + sigma_n*eye(length(Lambda(:,1)));
    C_inv = C^(-1);
    
    % Calculating the marginal likelihood multiplied by -1
    n = length(y);
    f = (0.5)*y'*C_inv*y + 0.5*log(norm(C)) + (n/2)*log(2*pi); 

end

