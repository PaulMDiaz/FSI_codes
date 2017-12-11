clear all
close all
clc

load('results_cylinder.mat')

% p_diff, cd, cl, drag, lift
p_diff_c = results(:,1); 
C_D_c = results(:,2); 
C_L_c = results(:,3); 
drag_c = results(:,4); 
lift_c = results(:,5); 
t_vec_c = linspace(0,6,length(lift_c))'; 




load('/Users/felixnewberry/Google Drive/1_Research/8_FSI/FSI_codes/FSI_cylinder/results_64.mat')

% ux, uy, drag lift

ux = results(:,1); 
uy = results(:,2); 
drag_fsi = results(:,3); 
lift_fsi = results(:,4); 

ux_true_l = 2.13e-5; 
ux_true_u = 2.27e-5; 

uy_true_u = 8.33e-4; 
uy_true_l = 8.16e-4; 

F_d_true_l = 14.2263; 
F_d_true_u = 14.38; 

F_l_true_l = 0.7517; 
F_l_true_u = 0.76487; 

% Calculate lift and drag coefficients. multiply by 2/(u_mean^2*2*0.05) 
u_mean = 0.2; 

C_D_fsi = 2/(u_mean^2*2*0.05)*drag_fsi; 
C_L_fsi = 2/(u_mean^2*2*0.05)*lift_fsi; 

t_vec_fsi = linspace(0,6,length(lift_fsi))'; 


% figure
% hold on
% h2 = plot(t_vec_c,drag_c,'-b', 'LineWidth', 2);
% h3 = plot(t_vec_fsi,drag_fsi,'-g', 'LineWidth', 2);
% ylim([0.01, 0.016]);
% 
% 
% figure
% hold on
% h2 = plot(t_vec_c,lift_c,'-b', 'LineWidth', 2);
% h3 = plot(t_vec_fsi,lift_fsi,'-g', 'LineWidth', 2);
% ylim([0.00, 14e-4]);

figure
hold on
h1 = plot(t_vec_fsi,F_d_true_l*ones(length(t_vec_fsi),1),'-r', 'LineWidth', 2);
h1 = plot(t_vec_fsi,F_d_true_u*ones(length(t_vec_fsi),1),'-r', 'LineWidth', 2);
h2 = plot(t_vec_c,C_D_c,'-b', 'LineWidth', 2);
h3 = plot(t_vec_fsi,C_D_fsi,'-g', 'LineWidth', 2);
legend([h1,h2,h3], {'Turek', 'Cylinder', 'FSI'},'interpreter', ...
        'latex', 'fontsize', 16);
hold off
xlabel('time', 'interpreter', 'latex', 'fontsize', 20);
ylabel('$C_D$', 'interpreter', 'latex', 'fontsize', 20);
ylim([0, 15])

figure
hold on
h1 = plot(t_vec_fsi,F_l_true_l*ones(length(t_vec_fsi),1),'-r', 'LineWidth', 2);
h1 = plot(t_vec_fsi,F_l_true_u*ones(length(t_vec_fsi),1),'-r', 'LineWidth', 2);
h2 = plot(t_vec_c,C_L_c,'-b', 'LineWidth', 2);
h3 = plot(t_vec_fsi,C_L_fsi,'-g', 'LineWidth', 2);
legend([h1,h2,h3], {'Turek', 'Cylinder', 'FSI'},'interpreter', ...
        'latex', 'fontsize', 16);
hold off
xlabel('time', 'interpreter', 'latex', 'fontsize', 20);
ylabel('$C_L$', 'interpreter', 'latex', 'fontsize', 20);
ylim([0.00, 0.8]);

figure 
hold on  
h1 = plot(t_vec_fsi,ux_true_u*ones(length(t_vec_fsi),1),'-r', 'LineWidth', 2);
h1 = plot(t_vec_fsi,ux_true_l*ones(length(t_vec_fsi),1),'-r', 'LineWidth', 2);
h2 = plot(t_vec_fsi,ux,'-b', 'LineWidth', 2);
legend([h1,h2], {'Turek', 'FSI'},'interpreter', ...
        'latex', 'fontsize', 16);
xlabel('time', 'interpreter', 'latex', 'fontsize', 20);
ylabel('$u_x$', 'interpreter', 'latex', 'fontsize', 20);
ylim([2e-5, 2.4e-5]);

figure 
hold on
h1 = plot(t_vec_fsi,uy_true_u*ones(length(t_vec_fsi),1),'-r', 'LineWidth', 2);
h1 = plot(t_vec_fsi,uy_true_l*ones(length(t_vec_fsi),1),'-r', 'LineWidth', 2);
h2 = plot(t_vec_fsi,uy,'-b', 'LineWidth', 2);
legend([h1,h2], {'Turek', 'FSI'},'interpreter', ...
        'latex', 'fontsize', 16);
xlabel('time', 'interpreter', 'latex', 'fontsize', 20);
ylabel('$u_y$', 'interpreter', 'latex', 'fontsize', 20);
ylim([-4e-4, 10e-4]);
