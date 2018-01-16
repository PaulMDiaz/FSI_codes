% Structure validation 
clear all
close all
clc

format long 

% load('results_7906_02.mat')
% load('results_7906_01.mat')
% load('results_7906_005.mat')

% load('results_29850_02.mat')
% load('results_29850_01.mat')
% load('results_29850_005.mat')

load('results_105134_02.mat') % treat special. 

% load('results_105134_01.mat')
% load('results_105134_005.mat')
% dt = 0.02
% t_1 = 400; t_2 = 450; t_3 = 500; t_0 = 50; 

% dt = 0.01
% t_1 = 800; t_2 = 900; t_3 = 1000; t_0 = 100;

% dt = 0.005
% t_1 = 1600; t_2 = 1800; t_3 = 2000; t_0 = 200;

t = results(1:end-1,1); 
u_x = results(1:end-1,2); 
u_y = results(1:end-1, 3); 


% Find average max and min in range from t = 4 to t = 4.8. changes for time
% step. 






x_max = max(u_x(t_1:t_2));

[x_min, min_loc_1] = min(u_x(t_1:t_2)); 
[x_min_1, min_loc_2] = min(u_x(t_2:t_3)); 
[x_min_2, min_loc_3] = min(u_x(1:t_0)); 

f_x = 10/(t(t_2+min_loc_2) - t(min_loc_3)); 

y_max = max(u_y(t_1:t_2));
[y_min, min_loc_1] = min(u_y(t_1:t_2)); 
[y_min_1, min_loc_2] = min(u_y(t_2:t_3)); 
[y_min_2, min_loc_3] = min(u_y(1:t_0)); 

f_y = 10/(t(t_2+min_loc_2) - t(min_loc_3)); 

x_var = (x_max - x_min)/2;
x_mean = (x_max + x_min)/2;
y_var = (y_max - y_min)/2;
y_mean = (y_max + y_min)/2;


[x_mean, x_var, f_x, y_mean, y_var, f_y]

% % Find mean, amplitude and frequency of u_x and u_y 
% Fs = 100;   % Sampling freq or 50, 100, 200
% T = 10;     % Sampling period 
% L = 1000    % Length of signal
% 
% n = 2^nextpow2(L);
% X = fft(u_y'); 
% 
% P2 = abs(X/n);
% P1 = P2(:,1:n/2+1);
% P1(:,2:end-1) = 2*P1(:,2:end-1);
[peaks, locations] = findpeaks(u_x); %all peaks and locations

figure
subplot(2,2,1)
plot(t, u_x,'LineWidth', 2)
xlabel('time', 'interpreter', 'latex', 'fontsize', 20);
ylabel('displacement x', 'interpreter', 'latex', 'fontsize', 20);
xlim([0, 10])
ylim([-0.03, 0.005])

subplot(2,2,2)
plot(t, u_y,'LineWidth', 2)
xlabel('time', 'interpreter', 'latex', 'fontsize', 20);
ylabel('displacement y', 'interpreter', 'latex', 'fontsize', 20);
xlim([0, 10])
ylim([-0.14, 0.02])

subplot(2,2,3)
plot(t, u_x,'LineWidth', 2)
xlabel('time', 'interpreter', 'latex', 'fontsize', 20);
ylabel('displacement x', 'interpreter', 'latex', 'fontsize', 20);
xlim([8, 10])
ylim([-0.03, 0.005])

subplot(2,2,4)
plot(t, u_y,'LineWidth', 2)
xlabel('time', 'interpreter', 'latex', 'fontsize', 20);
ylabel('displacement y', 'interpreter', 'latex', 'fontsize', 20);
xlim([8, 10])
ylim([-0.14, 0.02])
