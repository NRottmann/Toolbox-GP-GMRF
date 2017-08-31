function [obj] = random(num_all,num_cov,D,m,range,optimthings)
% This function creates random varibles as starting points for the
% minimization process. Call due to hypParam_opt.m

obj = zeros(num_all,1);
if optimthings == 0             % Defining Random start variables for different settings
    obj(1:num_cov) = rand(num_cov,1)*10;
    for ll = 1:1:D
        for jj=1:1:m
            obj(num_cov+(ll-1)*m+jj) = range(1,ll) + rand(1,1) * (range(2,ll) - range(1,ll));
        end
    end
elseif optimthings == 1
    obj = rand(num_cov,1)*range;
elseif optimthings == 2
    for ll = 1:1:D
        for jj=1:1:m
            obj((ll-1)*m+jj) = range(1,ll) + rand(1,1) * (range(2,ll) - range(1,ll));
        end
    end
else
    error('Wrong input for optimthings!')
end


end

