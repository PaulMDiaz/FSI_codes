clear all
close all
clc

format long 

load('nodal_disp.mat')
load('disp_coordinates2.mat')
% load('disp_coordinates3.mat')

coords_length = 410; %814; % 410 %294
interfaceStructureDispCoords = interfaceStructureDispCoords2; 

% nodal_disp
% 294 values

artifical_disp = zeros(coords_length,1); 

% try just looking at x coordinate. 

y_vec = interfaceStructureDispCoords(:,2); 

x_vec = interfaceStructureDispCoords(:,1); 
x_vec = x_vec-min(x_vec); 

% parabola to describe curve of structure
x_start = 0.24898979;
x_end = 0.6;

x_length = x_end - x_start; 

height = 0.005; 
disp_parabola = height/0.35^2*x_vec.^2; 

% I miss a little bit of length... 

% 0.35 is max length. 0.75 yields max of 0.091875



% sin curve fot time so that deflection oscillates. 

t = 0:0.001:8;
% t = 0:0.0005:8;

% artifical_disp = zeros(294,length(t)); 

% Sin curve in time
T = 1; 
disp_t = sin((2*pi)/(T)*t); 


% Sin^2 curve in time
% T = 1; 
% disp_t = sin((2*pi)/(T)*t).^2; 

% Step function. have displacement parabola apply at t = 5. 
T_step = 2; % time at which to apply step
% disp_t = ones(1,length(t)); 
% disp_t(1:T_step) = 0; 


figure
% plot(x_vec, disp_parabola)
plot(t, disp_t)
% xlim([0,2])

% Only adjust y coordiante, not x
disp_parabola(1:2:end) = 0; 

artifical_disp = disp_parabola*disp_t; 

% I think area is conserved
% need to save artificial_disp for use in python. 

% save('./artifical_disp.mat', 'artifical_disp');
% save('./artifical_disp_med2_h_005_stept2', 'artifical_disp','-ascii');

save('./artifical_disp_med2_h_005_parabola', 'disp_parabola','-ascii');
save('./artifical_disp_med2_h_005_time', 'disp_t','-ascii');

figure
plot(x_vec, y_vec)

% figure
% ylim([-0.1 0.1])
% for i = 1:5:length(t)
%     ylim([-0.01 0.01])
%     plot(x_vec(2:2:end), artifical_disp(2:2:end,i))
%     ylim([-0.01 0.01])
%     pause
%     hold off
% end
