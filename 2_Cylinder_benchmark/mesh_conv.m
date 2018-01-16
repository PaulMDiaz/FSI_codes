clear all
close all
clc

% Want to converge to C_D = 5.57 to 5.59, C_L 0.0104 to 0.0110, p_diff
% 0.1172 to 0.1176

% 32 has dofs_V 9888 and dofs_Q 1316
% 64 has dofs_V 37964 and dofs_Q 4908
% 128 is 147108, 18710. 
% 256 585792, 73868

% load('results_64_N.mat')
% load('results_cylinder.mat')
% load('results_cylinder_64_2N.mat')
load('results_cylinder_32.mat')

p_diff_32 = results(1:end-1,1); 
C_D_32 = results(1:end-1,2); 
C_L_32 = results(1:end-1,3); 
t_vec_32 = linspace(0,8,length(C_L_32))'; 

% load('results_64_N.mat')
load('results_cylinder_64.mat')
% load('results_cylinder_128.mat')
% this is really 64 with 2*N on cylinder. 

p_diff_64 = results(1:end-1,1); 
C_D_64 = results(1:end-1,2); 
C_L_64 = results(1:end-1,3); 
t_vec_64 = linspace(0,8,length(C_L_64))'; 

drag_64 = results(1:end-1,4); 
lift_64 = results(1:end-1,5); 

% load('results_cylinder_128.mat')
load('results_cylinder_96.mat')
p_diff_96 = results(1:end-1,1); 
C_D_96 = results(1:end-1,2); 
C_L_96 = results(1:end-1,3); 
t_vec_96 = linspace(0,8,length(C_L_96))'; 

load('results_cylinder_128.mat')
p_diff_128 = results(1:end-1,1); 
C_D_128 = results(1:end-1,2); 
C_L_128 = results(1:end-1,3); 
t_vec_128 = linspace(0,8,length(C_L_128))'; 
% 
% load('results_256.mat')
% p_diff_256 = results(:,1); 
% C_D_256 = results(:,2); 
% C_L_256 = -results(:,3); 
% t_vec_256 = linspace(0,8,length(C_L_256))'; 

CD_l = 5.57; 
CD_u = 5.59; 
CL_l = 0.0104;
CL_u = 0.0110; 
p_diff_l = 0.1172;
p_diff_u = 0.1176; 

figure
hold on
h1 = plot(t_vec_32, ones(1,length(t_vec_32))*CD_u,'-r', 'LineWidth', 2);
h1 = plot(t_vec_32, ones(1,length(t_vec_32))*CD_l,'-r', 'LineWidth', 2);
h2 = plot(t_vec_32,C_D_32,'-b', 'LineWidth', 2);
h3 = plot(t_vec_64,C_D_64,'-g', 'LineWidth', 2);
h4 = plot(t_vec_96,C_D_96,'-m', 'LineWidth', 2);
h5 = plot(t_vec_128,C_D_128,'-c', 'LineWidth', 2);
legend([h1,h2,h3,h4, h5], {'True', '32','64','96','128'},'interpreter', ...
        'latex', 'fontsize', 16);
    
hold off
xlabel('time', 'interpreter', 'latex', 'fontsize', 20);
ylabel('$C_D$', 'interpreter', 'latex', 'fontsize', 20);
ylim([5.5, 5.7])
% ylim([5, 13])

figure
hold on
h1 = plot(t_vec_32, ones(1,length(t_vec_32))*CL_u,'-r', 'LineWidth', 2);
h1 = plot(t_vec_32, ones(1,length(t_vec_32))*CL_l,'-r', 'LineWidth', 2);
h2 = plot(t_vec_32,C_L_32,'-b', 'LineWidth', 2);
h3 = plot(t_vec_64,C_L_64,'-g', 'LineWidth', 2);
h4 = plot(t_vec_96,C_L_96,'-m', 'LineWidth', 2);
h5 = plot(t_vec_128,C_L_128,'-c', 'LineWidth', 2);
legend([h1,h2,h3,h4, h5], {'True', '32','64','96','128'},'interpreter', ...
        'latex', 'fontsize', 16);
hold off
xlabel('time', 'interpreter', 'latex', 'fontsize', 20);
ylabel('$C_L$', 'interpreter', 'latex', 'fontsize', 20);
ylim([0, 0.04])
% ylim([0, 0.1])

figure
hold on
h1 = plot(t_vec_32, ones(1,length(t_vec_32))*p_diff_u,'-r', 'LineWidth', 2);
h1 = plot(t_vec_32, ones(1,length(t_vec_32))*p_diff_l,'-r', 'LineWidth', 2);
h2 = plot(t_vec_32,p_diff_32,'-b', 'LineWidth', 2);
h3 = plot(t_vec_64,p_diff_64,'-g', 'LineWidth', 2);
h4 = plot(t_vec_96,p_diff_96,'-m', 'LineWidth', 2);
h5 = plot(t_vec_128,p_diff_128,'-c', 'LineWidth', 2);
legend([h1,h2,h3,h4, h5], {'True', '32','64','96','128'},'interpreter', ...
        'latex', 'fontsize', 16);
hold off
xlabel('time', 'interpreter', 'latex', 'fontsize', 20)
ylabel('$p_{diff}$', 'interpreter', 'latex', 'fontsize', 20)
ylim([0.115, 0.125])

% Make table of results
% 32, 64, 96, 128 corresponds to # unknowns: 

% % my results
% L_me = [C_L_32(end), C_L_64(end), C_L_96(end), C_L_128(end)]; 
% D_me = [C_D_32(end), C_D_64(end), C_D_96(end), C_D_128(end)]; 
% P_me = [p_diff_32(end), p_diff_64(end), p_diff_96(end), p_diff_128(end)];
% 
% ndof_me = [5240+707, 18416+2406, 38842+4988, 69136+8850]; 
% 
% ndof_12 = [32592, 26970, 22212]; 
% D_12 = [5.569, 5.5125, 5.6026]; 
% L_12 = [0.0132, 0.0056, -0.0031]; 
% P_12 = [0.1155, 0.1154, 0.1167]; 
% 
% ndof_10 = [167232]; 
% D_12 = [5.569, 5.5125, 5.6026]; 
% L_12 = [0.0132, 0.0056, -0.0031]; 
% P_12 = [0.1155, 0.1154, 0.1167]; 
% % Plot against turek results
% 
% figure
% hold on
% h1 = plot(ndof_me, D_me,'-k*', 'LineWidth', 2);
% h12 = plot(ndof_12,D_12,'-bo', 'LineWidth', 2);
% 
% % legend([h1,h2,h3,h4, h5], {'True', '32','64','96','128'},'interpreter', ...
% %         'latex', 'fontsize', 16);
% hold off
% xlabel('ndofs', 'interpreter', 'latex', 'fontsize', 20)
% ylabel('$C_D$', 'interpreter', 'latex', 'fontsize', 20)
% % ylim([0.115, 0.125])
% 
% figure
% hold on
% h1 = plot(ndof_me, L_me,'-k*', 'LineWidth', 2);
% h12 = plot(ndof_12,L_12,'-bo', 'LineWidth', 2);
% 
% % legend([h1,h2,h3,h4, h5], {'True', '32','64','96','128'},'interpreter', ...
% %         'latex', 'fontsize', 16);
% hold off
% xlabel('ndofs', 'interpreter', 'latex', 'fontsize', 20)
% ylabel('$C_L$', 'interpreter', 'latex', 'fontsize', 20)
% % ylim([0.115, 0.125])
% 
% figure
% hold on
% h1 = plot(ndof_me, P_me,'-k*', 'LineWidth', 2);
% h12 = plot(ndof_12,P_12,'-bo', 'LineWidth', 2);
% 
% % legend([h1,h2,h3,h4, h5], {'True', '32','64','96','128'},'interpreter', ...
% %         'latex', 'fontsize', 16);
% hold off
% xlabel('ndofs', 'interpreter', 'latex', 'fontsize', 20)
% ylabel('$P_{diff}$', 'interpreter', 'latex', 'fontsize', 20)
% % ylim([0.115, 0.125])
