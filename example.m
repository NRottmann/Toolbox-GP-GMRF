% Script to evaluate the GP and GMRF for a stationary, created field.
clc
close all
clear all

% include folder others
Full = mfilename('fullpath');
Name = mfilename;
base = Full(1:end-length(Name));
addpath([base '\others']);


%%
% get the data
M = open('field.mat');

x_num = length(M.LonScale);                             % Number of steps in x-direction
y_num = length(M.LatScale);                             % Number of steps in x-direction

x_range = M.LonScale(x_num) - M.LonScale(1);    % range in x-direction [°]   
y_range = M.LatScale(y_num) - M.LatScale(1);    % range in y-direction [°]

Field = M.TempField(:,:,1);                     % Temperature Field

%%
%%%%%%%%%%%%%%%%%
% Generating the estimation points (s) und generating points (p)
%%%%%%%%%%%%%%%%%
% real data mesh points
num_s = 100;
num_p = 5;

[s,~] = MeshGen([num_s; num_s],[x_range; y_range]);
[p,~] = MeshGen([num_p; num_p],[x_range; y_range]);

%%
%%%%%%%%%%%%%%%%%%%%
% Making measurements and calculating new field estimations
%%%%%%%%%%%%%%%%%%%%
nn = 30;                                % Number of measuring points
x = zeros(2,nn);
for jj = 1:1:nn
    x_measure = rand(1,1) * x_range;      % creating random measure positions at time t = i
    y_measure = rand(1,1) * y_range;
    while isnan(get_f_bil([x_measure; y_measure],Field,[x_range; y_range]))
        x_measure = rand(1,1) * x_range;     
        y_measure = rand(1,1) * y_range;
    end
    x(:,jj) = [x_measure; y_measure];
end

% Get the values at the measuring positions X_measure
[y] = get_f_bil(x,Field,[x_range; y_range]);

%%
% GMRF using the created toolbox
p = [2; 30; 0; x_range; 0; y_range];
[mu,Sigma] = GMRF_Spatial(p,x,s,y,'noise',0.1,'improve',true,'method',1,'maxTime',1);


%%
% Plotting the field
Field_mean = zeros(num_s,num_s);
Field_variance = zeros(num_s,num_s);

% check for nan's
[y_s] = get_f_NN(s,Field,[x_range; y_range]);
II = isnan(y_s);
mu(II) = nan;
Sigma(II) = nan;

for ii=1:1:num_s
    Field_mean(:,ii) = mu((1+num_s*(ii-1)):(num_s*ii));
    Field_variance(:,ii) = Sigma((1+num_s*(ii-1)):(num_s*ii));
end

max_value = max(max(Field)) + 1;   % Defining the maximum values for the colorbar
min_value = min(min(Field)) - 1;

figure(1)
subplot(1,2,1)
hold on;
plot_field(Field,20,[x_range; y_range],max_value,min_value)
title('Original Field')

subplot(1,2,2)
hold on;
plot_field(Field_mean,20,[x_range; y_range],max_value,min_value)
hold on;
plot(x(1,:),x(2,:),'k*')
hold off;