function [f,g] = logLikelihood_GMRF(p_ord,x,y,sigma_n,CovFunc,CovParam,PreFunc,PreParam,D,m,optimthings)
    % Loglikelihood Function for GP with built-in GMRF
    % This is the Loglikelihood function which has to be minimized for
    % finding optimal parameter for a GP with built-in GMRF. Call due to
    % hypParam_opt.m


    % define functions to call them
    Pre = str2func(PreFunc);
    Cov = str2func(CovFunc);

    % generating points have a predefined ordering
    for c=1:1:D
        p(c,:) = p_ord((1+(m*(c-1))):(c*m));
    end
    
    % Calculation of the COvariance Matrix and the precision matrix
    [Lambda] = Cov(x,p,'CovParam',CovParam);

    % Calculating the Covariance Matrix for y
    if sigma_n > 0.001
        [Q] = Pre(p,'PreParam',PreParam);
        C = Lambda*Q^(-1)*Lambda' + sigma_n*eye(length(Lambda(:,1)));
    else
        C = Lambda * Lambda' + 0.1 * eye(length(Lambda(:,1)));
    end
    C_inv = inv(C);
    
    % Calculating the marginal likelihood multiplied by -1
    n = length(y);
    f = (0.5)*y'*C_inv*y + 0.5*log(norm(C)) + (n/2)*log(2*pi); 

    % Calculation of the derivation of the Matrix K with regard to the
    % hyperparameter using (f(x_0 + h) - f(x)) / h
    alpha = C_inv * y;      
    h = 0.001;          % step size for covariance function
    
    if optimthings == 0
        num_p = length(p_ord);
        num_cov = length(CovParam);
        num_all = num_p + num_cov;
        g = zeros(num_all,1);
        for i=1:1:num_cov
            CovParam_h = CovParam;                  % Calculating a change in Hyperparameter
            CovParam_h(i) = CovParam_h(i) + h;
            [Lambda_h] = Cov(x,p,'CovParam',CovParam_h);
            if sigma_n > 0.001
                C_h = Lambda_h*Q^(-1)*Lambda_h' + sigma_n*eye(length(Lambda_h(:,1))); % Calculating the partial derivative
            else
                C_h = Lambda_h * Lambda_h' + 0.1 * eye(length(Lambda_h(:,1)));
            end 
            C_d = ((C_h - C)/h);
            g(i) = (0.5*trace((C_inv - alpha*alpha')*C_d));
        end
        for i=1:1:num_p
            p_ord_h = p_ord;                      % Calculating a change in Hyperparameter
            p_ord_h(i) = p_ord_h(i) + h;
            for c=1:1:D
                p_h(c,:) = p_ord_h((1+(m*(c-1))):(c*m));
            end
            [Lambda_h] = Cov(x,p_h,'CovParam',CovParam);
            if sigma_n > 0.001
                [Q_h] = Pre(p_h,'PreParam',PreParam);
                C_h = Lambda_h*Q_h^(-1)*Lambda_h' + sigma_n*eye(length(Lambda_h(:,1)));      % Calculating the partial derivative
            else
                C_h = Lambda_h * Lambda_h' + 0.1 * eye(length(Lambda_h(:,1)));
            end
            C_d = ((C_h - C)/h);
            g(num_cov+i) = (0.5*trace((C_inv - alpha*alpha')*C_d));
        end
    elseif optimthings == 1
        num_cov = length(CovParam);
        g = zeros(num_cov,1);
        for i=1:1:num_cov
            CovParam_h = CovParam;                  % Calculating a change in Hyperparameter
            CovParam_h(i) = CovParam_h(i) + h;
            [Lambda_h] = Cov(x,p,'CovParam',CovParam_h);
            if sigma_n > 0.001
                C_h = Lambda_h*Q^(-1)*Lambda_h' + sigma_n*eye(length(Lambda_h(:,1)));       % Calculating the partial derivative
            else
                C_h = Lambda_h * Lambda_h' + 0.1 * eye(length(Lambda_h(:,1)));
            end
            C_d = ((C_h - C)/h);
            g(i) = (0.5*trace((C_inv - alpha*alpha')*C_d));
        end
    elseif optimthings == 2
        num_p = length(p_ord);
        g = zeros(num_p,1);
        for i=1:1:num_p
            p_ord_h = p_ord;                      % Calculating a change in Hyperparameter
            p_ord_h(i) = p_ord_h(i) + h;
            for c=1:1:D
                p_h(c,:) = p_ord_h((1+(m*(c-1))):(c*m));
            end
            [Lambda_h] = Cov(x,p_h,'CovParam',CovParam);
            if sigma_n > 0.001
                [Q_h] = Pre(p_h,'PreParam',PreParam);
                C_h = Lambda_h*Q_h^(-1)*Lambda_h' + sigma_n*eye(length(Lambda_h(:,1)));      % Calculating the partial derivative
            else
                C_h = Lambda_h * Lambda_h' + 0.1 * eye(length(Lambda_h(:,1)));
            end
            C_d = ((C_h - C)/h);
            g(i) = (0.5*trace((C_inv - alpha*alpha')*C_d));
        end
    end
end

