% clear all
% close all
% clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Convergence to Turek FSI benchmark %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Need to fix code before this is relevant. 

% % load('results_1.mat')
% load('cylinder_fsi_course_full.mat')
% 
% 
% 
% dt = 0.001; 
% u1_64 = results(1:end-1,1); 
% u2_64 = results(1:end-1,2); 
% drag_64 = results(1:end-1,3); 
% lift_64 = results(1:end-1,4); 
% t_vec_64 = [0:dt:dt*(length(lift_64)-1)];

load('cylinder_fsi_course.mat')

dt = 0.001; 
u1_1 = results(1:end-1,1); 
u2_1 = results(1:end-1,2); 
drag_1 = results(1:end-1,3); 
lift_1 = results(1:end-1,4); 
t_vec_1 = [0:dt:dt*(length(lift_1)-1)];


% load('cylinder_bar_fluid_fine.mat')
% load('cylinder_bar_fluid_med.mat')
% 
% 
% drag_no_traction = results(1:end-1,1); 
% lift_no_traction = results(1:end-1,2); 
% t_vec_no_traction = [0:dt:dt*(length(lift_no_traction)-1)];

% figure
% plot(t_vec_64, lift_64)
% 
% figure
% plot(t_vec_64, drag_64)


% load('results_64_structure_001.mat')
% u1 = results(:,2); 
% u2 = results(:,3); 
% dt = 0.001; 
% t_struct = 0:dt:dt*(length(u1)-1); 
% 
% load('results_64_structure_FSI.mat')
% u1_fsi = results(:,1); 
% u2_fsi = results(:,2); 
% dt_fsi = 0.001
% t_struct_fsi = 0:dt_fsi:dt_fsi*(length(u1_fsi)-1); 
% 
% figure
% hold on
% plot(t_struct, u1)
% plot(t_struct_fsi, u1_fsi)
% 
% figure
% hold on
% plot(t_struct, u2)
% plot(t_struct_fsi, u2_fsi)

% t_vec_64 = linspace(0,10,length(u1_64))'; 
% 
u1 = 2.27e-5;
u2 = 8.209e-4; 
drag_fsi = 14.295; 
lift_fsi = 0.7638; 

drag_cfd = 14.29; 
lift_cfd = 1.119; 
% 
figure
hold on
h1 = plot(t_vec_1,u1*ones(1,length(t_vec_1)),'-r', 'LineWidth', 2);
% h2 = plot(t_vec_64,u1_64,'-b', 'LineWidth', 2);
h3 = plot(t_vec_1,u1_1,'-k', 'LineWidth', 2);
% legend([h1,h2], {'True', 'FSI course'},'interpreter', ...
%         'latex', 'fontsize', 16);
xlabel('time', 'interpreter', 'latex', 'fontsize', 20);
ylabel('$u_x$', 'interpreter', 'latex', 'fontsize', 20);

ylim([1.4e-5, 2.4e-5])
% % ylim([0, 0.1])
% 

figure
hold on
h1 = plot(t_vec_1,u2*ones(1,length(t_vec_1)),'-r', 'LineWidth', 2);
% h2 = plot(t_vec_64,u2_64,'-b', 'LineWidth', 2);
h3 = plot(t_vec_1,u2_1,'-k', 'LineWidth', 2);
% legend([h1,h2], {'True', 'FSI course'},'interpreter', ...
%         'latex', 'fontsize', 16);
xlabel('time', 'interpreter', 'latex', 'fontsize', 20);
ylabel('$u_y$', 'interpreter', 'latex', 'fontsize', 20);
ylim([0, 10e-4])

figure
hold on
h1 = plot(t_vec_1,drag_fsi*ones(1,length(t_vec_1)),'-r', 'LineWidth', 2);
% h2 = plot(t_vec_64,drag_64,'-b', 'LineWidth', 2);
h5 = plot(t_vec_1,drag_1,'-k', 'LineWidth', 2);
h3 = plot(t_vec_1,drag_cfd*ones(1,length(t_vec_1)),'--r', 'LineWidth', 2);
% h4 = plot(t_vec_no_traction,drag_no_traction,'-g', 'LineWidth', 2);

% legend([h1,h2, h3, h4], {'True FSI', 'FSI_course','True CFD','CFD_med'},'interpreter', ...
%         'latex', 'fontsize', 16);
xlabel('time', 'interpreter', 'latex', 'fontsize', 20);
ylabel('drag', 'interpreter', 'latex', 'fontsize', 20);
ylim([0, 15])

figure
hold on
h1 = plot(t_vec_1,lift_fsi*ones(1,length(t_vec_1)),'-r', 'LineWidth', 2);
% h2 = plot(t_vec_64,lift_64,'-b', 'LineWidth', 2);
h5 = plot(t_vec_1,lift_1,'-k', 'LineWidth', 2);
h3 = plot(t_vec_1,lift_cfd*ones(1,length(t_vec_1)),'--r', 'LineWidth', 2);
% h4 = plot(t_vec_no_traction,lift_no_traction,'-g', 'LineWidth', 2);

% legend([h1,h2, h3, h4], {'True FSI', 'FSI course','True CFD','CFD_med'},'interpreter', ...
%         'latex', 'fontsize', 16);
xlabel('time', 'interpreter', 'latex', 'fontsize', 20);
ylabel('lift', 'interpreter', 'latex', 'fontsize', 20);
ylim([-01, 1.4])
