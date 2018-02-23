clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Convergence to Turek FSI benchmark %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Need to fix code before this is relevant. 


% load('results_32_fluid.mat')

% dt = 0.001; 
% u1_64 = results(1:end-1,1); 
% u2_64 = results(1:end-1,2); 
% drag_64 = results(1:end-1,3); 
% lift_64 = results(1:end-1,4); 
% t_vec_64 = [0:dt:dt*(length(lift_64)-1)];

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
% u1 = 2.2630e-5;
% u2 = 8.2935e-4; 
% drag = 14.2970; 
% lift = 0.76687; 
% 
% figure
% hold on
% h1 = plot(t_vec_64,u1*ones(1,length(t_vec_64)),'-r', 'LineWidth', 2);
% h2 = plot(t_vec_64,u1_64,'-b', 'LineWidth', 2);
% legend([h1,h2], {'True', '64'},'interpreter', ...
%         'latex', 'fontsize', 16);
% xlabel('time', 'interpreter', 'latex', 'fontsize', 20);
% ylabel('$u_x$', 'interpreter', 'latex', 'fontsize', 20);

% % ylim([0, 0.1])
% % % ylim([0, 0.1])
% % 
% 
% figure
% hold on
% h1 = plot(t_vec_64,u2*ones(1,length(t_vec_64)),'-r', 'LineWidth', 2);
% h2 = plot(t_vec_64,u2_64,'-b', 'LineWidth', 2);
% legend([h1,h2], {'True', '64'},'interpreter', ...
%         'latex', 'fontsize', 16);
% xlabel('time', 'interpreter', 'latex', 'fontsize', 20);
% ylabel('$u_y$', 'interpreter', 'latex', 'fontsize', 20);
% 
% figure
% hold on
% h1 = plot(t_vec_64,drag*ones(1,length(t_vec_64)),'-r', 'LineWidth', 2);
% h2 = plot(t_vec_64,drag_64,'-b', 'LineWidth', 2);
% legend([h1,h2], {'True', '64'},'interpreter', ...
%         'latex', 'fontsize', 16);
% xlabel('time', 'interpreter', 'latex', 'fontsize', 20);
% ylabel('drag', 'interpreter', 'latex', 'fontsize', 20);
% 
% figure
% hold on
% h1 = plot(t_vec_64,lift*ones(1,length(t_vec_64)),'-r', 'LineWidth', 2);
% h2 = plot(t_vec_64,lift_64,'-b', 'LineWidth', 2);
% legend([h1,h2], {'True', '64'},'interpreter', ...
%         'latex', 'fontsize', 16);
% xlabel('time', 'interpreter', 'latex', 'fontsize', 20);
% ylabel('lift', 'interpreter', 'latex', 'fontsize', 20);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Test Boundary Conditions: Traction %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('traction_tensor.mat')
% size is # points along FSI, measured quantities, t steps. 
% Columns are f_Traction_x,f_Traction_y,s_Traction_x,s_Traction_y
% where f is fluid traction force and s is structure
% traction_tensor = traction_tensor(:,:,1:1147);
t_print = 1;
t_length = length(traction_tensor(1,1,:));

% currently have 1147; 
% put in a pause and print every time step. 
figure
for i_traction = 1: t_length;
    t_print = i_traction; 
    
%     hold on 
    subplot(1,2,1)
    h1 = plot(traction_tensor(:,1,t_print), '-r', 'LineWidth', 2);
    hold on         
    h2 = plot(-traction_tensor(:,3,t_print), '--b', 'LineWidth', 2);
    hold off
%     legend([h1,h2], {'fluid','structure'},'interpreter', ...
%             'latex', 'fontsize', 16);
    
    xlabel('FSI nodal point', 'interpreter', 'latex', 'fontsize', 20)
    ylabel('$F_{x}$', 'interpreter', 'latex', 'fontsize', 20)
%     ylim([-100, 200]);
    
    subplot(1,2,2)
    h1 = plot(traction_tensor(:,2,t_print), '-r', 'LineWidth', 2);
    hold on 
    h2 = plot(-traction_tensor(:,4,t_print), '--b', 'LineWidth', 2);
    hold off
    legend([h1,h2], {'$T_{fluid}$','$T_{structure}$'},'interpreter', ...
            'latex', 'fontsize', 16);
    
    xlabel('FSI nodal point', 'interpreter', 'latex', 'fontsize', 20)
    ylabel('$F_{y}$', 'interpreter', 'latex', 'fontsize', 20)
%     ylim([-500, -200]); 
    title(['Time: ' num2str(t_print)])

    pause
end 

% hold off 
% ylim([0.11, 0.135])
% 
% figure
% subplot(1,2,2)
% hold on 
% h1 = plot(traction_tensor(:,2,t_print), '-r', 'LineWidth', 2);
% h2 = plot(-traction_tensor(:,4,t_print), '--b', 'LineWidth', 2);
% legend([h1,h2], {'fluid','structure'},'interpreter', ...
%         'latex', 'fontsize', 16);
% hold off
% xlabel('FSI nodal point', 'interpreter', 'latex', 'fontsize', 20)
% ylabel('$F_{y}$', 'interpreter', 'latex', 'fontsize', 20)
% % ylim([0.11, 0.135])

% load fluid, mesh and structure velocites. should all match. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Test Boundary Conditions: Fluid, Structure and Mesh velocities %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('u_m_FSI.mat')
load('v_s_FSI.mat')
load('u_f_FSI.mat')



% seperate x and y. 
u_m_x = u_m_FSI(1:2:end);
u_m_y = u_m_FSI(2:2:end);
v_s_x = v_s_FSI(1:2:end);
v_s_y = v_s_FSI(2:2:end);
u_f_x = u_f_FSI(1:2:end);
u_f_y = u_f_FSI(2:2:end);

%reverse top row for continuity of nodes
% This only works for the given mesh so not an excellent solution
% u_m_x(141:end) = u_m_x(end:-1:141); 
% u_m_y(141:end) = u_m_y(end:-1:141); 
% v_s_x(141:end) = v_s_x(end:-1:141); 
% v_s_y(141:end) = v_s_y(end:-1:141); 
% u_f_x(141:end) = u_f_x(end:-1:141); 
% u_f_y(141:end) = u_f_y(end:-1:141);
n_reverse = 141;
u_m_x(n_reverse:end) = u_m_x(end:-1:n_reverse); 
u_m_y(n_reverse:end) = u_m_y(end:-1:n_reverse); 
v_s_x(n_reverse:end) = v_s_x(end:-1:n_reverse); 
v_s_y(n_reverse:end) = v_s_y(end:-1:n_reverse); 
u_f_x(n_reverse:end) = u_f_x(end:-1:n_reverse); 
u_f_y(n_reverse:end) = u_f_y(end:-1:n_reverse);
% Nodes for this plot go bottom left to right, up the end then jump to top
% left to right. (when mesh in question is there... )

figure
subplot(1,2,1)
hold on 
plot(u_m_x, 'r', 'Linewidth', 2)
plot(u_f_x, '-.k', 'Linewidth', 2)
plot(v_s_x, '--b', 'Linewidth', 2)
xlabel('FSI nodes')
ylabel('x velocity')
hold off

subplot(1,2,2)
hold on 
h1 = plot(u_m_y, 'r', 'Linewidth', 2);
h2 = plot(u_f_y, '-.k', 'Linewidth', 2);
h3 = plot(v_s_y, '--b', 'Linewidth', 2);
xlabel('FSI nodes')
ylabel('y velocity')
hold off
legend([h1,h2, h3], {'$u_{mesh}$','$u_{fluid}$', '$u_{struct}$'},'interpreter', ...
            'latex', 'fontsize', 16);

