function [f,g] = BayOpt_objFun(x,s,y,num_all)
% Objective function for Bayesian Optimization
% This is the objective for the Bayesian Optimization which has to be
% minimized to find optimal parameters. Call with hypParam_opt.m

% parameter for bayesian optimization
alpha = 100;
h = 0.001;

% function value
[mu,Sigma] = GP_opt(x,s,y);
f = - mu - alpha * Sigma;
    
% differentiation of f numerically
g = zeros(1,num_all);
for ii=1:1:num_all
    s_h = s;                      
    s_h(ii) = s_h(ii) + h;        
    [mu_h,sigma_h] = GP_opt(x,s_h,y);
    f_h = - mu_h - alpha*sigma_h;
    g(ii) = ((f_h - f)/h);
end
end




