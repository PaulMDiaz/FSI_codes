clear all
close all
clc

load('nodal_disp.mat')
load('disp_coordinates.mat')

% interfaceStructureDispCoords

% nodal_disp
% 294 values

artifical_disp = zeros(294,1); 

% try just looking at x coordinate. 

y_vec = interfaceStructureDispCoords(:,2); 

x_vec = interfaceStructureDispCoords(:,1); 
x_vec = x_vec-min(x_vec); 

% parabola to describe curve of structure
height = 0.001; 
disp_parabola = height/0.35^2*x_vec.^2; 

% 0.35 is max length. 0.75 yields max of 0.091875


% sin curve fot time so that deflection oscillates. 

t = 0:0.001:8;
% artifical_disp = zeros(294,length(t)); 

% Sin curve in time
T = 32; 
disp_t = sin((2*pi)/(T)*t); 


% Sin^2 curve in time
T = 32; 
disp_t = sin((2*pi)/(T)*t).^2; 

% Step function. have displacement parabola apply at t = 5. 
T_step = 5; % time at which to apply step
disp_t = ones(1,length(t)); 
disp_t(1:5) = 0; 


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
save('./artifical_disp_step_000001', 'artifical_disp','-ascii');

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
