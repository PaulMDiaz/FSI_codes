% Structure validation 
clear all
close all
clc

format long 

% Load results

% load('results_7906_02.mat')
% load('results_7906_01.mat')
% load('results_7906_005.mat')

% load('results_29850_02.mat')
% load('results_29850_01.mat')
% load('results_29850_005.mat')

% load('results_105134_02.mat') 
% load('results_105134_01.mat')
load('results_105134_005.mat')

% Set dt at which data was generated:
% dt = 0.02;
% dt = 0.01;
dt = 0.005;

% Results contain time, u_x and u_y 
t = results(1:end-1,1); 
u_x = results(1:end-1,2); 
u_y = results(1:end-1, 3); 

% Turek et al. Caluclate mean and amplitude based on final period of
% oscillation. 

% Find max and min in final period. 
% final period can be taken as from 9 until end. 
u_max = max(u_x(9/dt:end));
u_min = min(u_x(9/dt:end));
y_max = max(u_y(9/dt:end));
y_min = min(u_y(9/dt:end));

% Require mean, amplitude and frequency. 

% mean in final period.  
mean_x = (u_max + u_min)/2;
mean_y = (y_max + y_min)/2;

% amplitude in final period. 
amp_x = (u_max - u_min)/2;
amp_y = (y_max - y_min)/2;

% frequency can do final preiod 1/T. 

% % Calculate upper and lower peaks (denoted by u and l) in x and y. MPD to prevent local max/min
MPD = 0.5/dt; 
[peaks_x_u, locations_x_u] = findpeaks(u_x,'MinPeakDistance',MPD); 
[peaks_x_l, locations_x_l] = findpeaks(-u_x,'MinPeakDistance',MPD); 
[peaks_y_u, locations_y_u] = findpeaks(u_y,'MinPeakDistance',MPD); 
[peaks_y_l, locations_y_l] = findpeaks(-u_y,'MinPeakDistance',MPD); 

freq_end_l = 1/(t(locations_x_l(end))-t(locations_x_l(end-1)));
freq_end_u = 1/(t(locations_x_u(end))-t(locations_x_u(end-1)));

% average of all periods 1/T
freq_all_lx = 10/(t(locations_x_l(end)) - t(locations_x_l(1))); 
freq_all_ly = 10/(t(locations_y_l(end)) - t(locations_y_l(1))); 
freq_all_u = 9/(t(locations_x_u(end)) - t(locations_x_u(1))); 

% Fourier Transform
fft(u_x);
Fs = 1/dt; % Sampling freq
x = u_x;
x = detrend(x,0); % take out regression line for more accurate fft
xdft = fft(x);
L = length(x); 

% freq = 0:Fs/length(x):Fs/2;
% I only need to search 1/2 of xdft for the max because x is real-valued
[Y,I] = max(abs(xdft(1:L/2+1)));
freq = 0:Fs/length(x):Fs/2;
freq_fft = freq(I);
% fprintf('Maximum occurs at %3.2f Hz\n.',freq(I))

% Plot freq
P2 = abs(xdft/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;
figure
plot(f,P1);
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')


% Print solution, mean, amplitude and freq. 
fprintf('x disp: %i +- %i [ % i ] \n', mean_x, amp_x, freq_all_lx);
fprintf('y disp: %i +- %i [ % i ] \n', mean_y, amp_y, freq_all_ly);

[mean_x, amp_x, mean_y, amp_y];
[freq_fft, freq_end_l, freq_end_u,freq_all_lx, freq_all_lx];

fprintf('%i', [freq_fft, freq_end_l, freq_end_u,freq_all_lx, freq_all_lx]);

% changing amplitude 
amp = (mean(u_x(locations_x_u) - u_x(locations_x_l(2:end))))/2;
x_bar = (mean(u_x(locations_x_u) + u_x(locations_x_l(2:end))))/2;

amp_change = (u_x(locations_x_u) - u_x(locations_x_l(2:end)))/2; 
x_bar_change = (u_x(locations_x_u) + u_x(locations_x_l(2:end)))/2; 
t_change = t(locations_x_u);




figure 
hold on
h1 = plot(t_change, x_bar_change, '-k', 'LineWidth', 2);
h2 = plot(t_change, x_bar_change + amp_change, '--k', 'LineWidth', 2);
h3 = plot(t_change, x_bar_change - amp_change, '--k', 'LineWidth', 2);
h4 = plot(t_change, 1e-3*(-14.305)*ones(length(t_change), 1), '-r', 'LineWidth', 1);
h5 = plot(t_change, 1e-3*(-14.305-14.305)*ones(length(t_change), 1), '-.r', 'LineWidth', 1);
h6 = plot(t_change, 1e-3*(-14.305+14.305)*ones(length(t_change), 1), '-.r', 'LineWidth', 1);
hold off
xlabel('time', 'interpreter', 'latex', 'fontsize', 20);
ylabel('displacement x', 'interpreter', 'latex', 'fontsize', 20);




% xlim([0, 10])
% ylim([-0.03, 0.005])
%all peaks and locations

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


figure
subplot(2,2,1)
hold on
plot(t, u_x,'LineWidth', 2)
plot(t(locations_x_u), u_x(locations_x_u), 'o')
plot(t(locations_x_l), u_x(locations_x_l), 'o')
hold off
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
% 
% x=u_x;
% Fs=1/dt;
% L=length(x);
% NFFT = 2^nextpow2(L);
% f = Fs/2*linspace(0,1,NFFT/2+1);
% 
% x=u_x;
% Fs = 1/dt; % Sampling freq
% L=length(x);
% NFFT = 2^nextpow2(L);
% f = Fs/2*linspace(0,1,NFFT/2+1);
% Y = fft(x,NFFT)/L;
% figure
% plot(f,2*abs(Y(1:NFFT/2+1))) 
% title('Single-Sided Amplitude Spectrum of y(t)')
% xlabel('Frequency (Hz)')
% ylabel('|Y(f)|')
% y_f = 2*abs(Y(1:NFFT/2+1)); 
% [a, I] = max(y_f(2:end));
% freq = f(I+1)
% 
% 
