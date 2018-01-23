clear all
close all
clc

% Want to converge to C_D = 5.57 to 5.59, C_L 0.0104 to 0.0110, p_diff
% 0.1172 to 0.1176


% 32 is V 7512 P 991
% 64 is V 25472 P 3288
% 128 is V 96008 P 12209
% 256 is V 369968 P 46662
n_dof = [7512+991, 25472+3288, 96008+12209, 369968+46662]';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('cylinder_chorin_32_01.mat')

p_diff_32 = results(1:end-1,1); 
C_D_32 = results(1:end-1,2); 
C_L_32 = results(1:end-1,3); 
t_vec_32 = linspace(0,15,length(C_L_32))'; 

load('cylinder_chorin_64_001.mat')
% 0.01 for N = 64 has nasty variations. 
% 0.005 is better. 

p_diff_64 = results(1:end-1,1); 
C_D_64 = results(1:end-1,2); 
C_L_64 = results(1:end-1,3); 
t_vec_64 = linspace(0,15,length(C_L_64))'; 

load('cylinder_chorin_128_0005.mat')

p_diff_128 = results(1:end-1,1); 
C_D_128 = results(1:end-1,2); 
C_L_128 = results(1:end-1,3); 
t_vec_128 = linspace(0,15,length(C_L_128))'; 

load('cylinder_chorin_128_0005.mat')

p_diff_256 = results(1:end-1,1); 
C_D_256 = results(1:end-1,2); 
C_L_256 = results(1:end-1,3); 
t_vec_256 = linspace(0,15,length(C_L_256))'; 

% Turek et al. Bounds. 
C_drag_l = 5.57; 
C_drag_u = 5.59; 
CD_average = (C_drag_l + C_drag_u)/2; 
C_lift_l = 0.0104;
C_lift_u = 0.0110; 
CL_average = (C_lift_l + C_lift_u)/2; 
p_diff_l = 0.1172;
p_diff_u = 0.1176; 
p_diff_average = (p_diff_l + p_diff_u)/2; 

figure
hold on
h1 = plot(t_vec_32, ones(1,length(t_vec_32))*C_drag_u,'-r', 'LineWidth', 2);
h1 = plot(t_vec_32, ones(1,length(t_vec_32))*C_drag_l,'-r', 'LineWidth', 2);
h2 = plot(t_vec_32,C_D_32,'-b', 'LineWidth', 2);
h3 = plot(t_vec_64,C_D_64,'-g', 'LineWidth', 2);
h4 = plot(t_vec_128,C_D_128,'-m', 'LineWidth', 2);
h5 = plot(t_vec_256,C_D_256,'-c', 'LineWidth', 2);
legend([h1,h2,h3,h4, h5], {'True', '32','64','128','256'},'interpreter', ...
        'latex', 'fontsize', 16);
    
hold off
xlabel('time', 'interpreter', 'latex', 'fontsize', 20);
ylabel('$C_D$', 'interpreter', 'latex', 'fontsize', 20);
ylim([5.5, 5.7])
% ylim([5, 13])

figure
hold on
h1 = plot(t_vec_32, ones(1,length(t_vec_32))*C_lift_u,'-r', 'LineWidth', 2);
h1 = plot(t_vec_32, ones(1,length(t_vec_32))*C_lift_l,'-r', 'LineWidth', 2);
h2 = plot(t_vec_32,C_L_32,'-b', 'LineWidth', 2);
h3 = plot(t_vec_64,C_L_64,'-g', 'LineWidth', 2);
h4 = plot(t_vec_128,C_L_128,'-m', 'LineWidth', 2);
h5 = plot(t_vec_256,C_L_256,'-c', 'LineWidth', 2);
legend([h1,h2,h3,h4, h5], {'True', '32','64','128','256'},'interpreter', ...
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
h4 = plot(t_vec_128,p_diff_128,'-m', 'LineWidth', 2);
h5 = plot(t_vec_256,p_diff_256,'-c', 'LineWidth', 2);
legend([h1,h2,h3,h4, h5], {'True', '32','64','128','256'},'interpreter', ...
        'latex', 'fontsize', 16);
hold off
xlabel('time', 'interpreter', 'latex', 'fontsize', 20)
ylabel('$p_{diff}$', 'interpreter', 'latex', 'fontsize', 20)
ylim([0.115, 0.125])

% Make table of results
% 32, 64,128, 256 corresponds to # unknowns  in n_dof. 

% my results
Lift_fenics = [C_L_32(end), C_L_64(end), C_L_128(end), C_L_256(end)]; 
Drag_fenics = [C_D_32(end), C_D_64(end), C_D_128(end), C_D_256(end)]; 
P_diff_fenics = [p_diff_32(end), p_diff_64(end), p_diff_128(end), p_diff_256(end)];

figure
xlim_vec = [8000, 1000000];
hold on
h2 = plot(n_dof, Drag_fenics,'-kx', 'LineWidth', 2);
h1 = plot(xlim_vec, ones(1,length(xlim_vec))*C_drag_u,'--r', 'LineWidth', 2);
h1 = plot(xlim_vec, ones(1,length(xlim_vec))*C_drag_l,'--r', 'LineWidth', 2);
hold off
legend([h1,h2], {'Turek et al.', 'Fenics_IPCS'},'interpreter', ...
        'latex', 'fontsize', 16, 'Location','Best');
set(gca, 'XScale', 'log')
xlim(xlim_vec);
xlabel('ndofs', 'interpreter', 'latex', 'fontsize', 20)
ylabel('$C_{Drag}$', 'interpreter', 'latex', 'fontsize', 20)

figure
hold on
h2 = plot(n_dof, Lift_fenics,'-k*', 'LineWidth', 2);
h1 = plot(xlim_vec, ones(1,length(xlim_vec))*C_lift_u,'--r', 'LineWidth', 2);
h1 = plot(xlim_vec, ones(1,length(xlim_vec))*C_lift_l,'--r', 'LineWidth', 2);
hold off
legend([h1,h2], {'Turek et al.', 'Fenics_IPCS'},'interpreter', ...
        'latex', 'fontsize', 16, 'Location','Best');
set(gca, 'XScale', 'log')
xlim(xlim_vec);
xlabel('ndofs', 'interpreter', 'latex', 'fontsize', 20)
ylabel('$C_{Lift}$', 'interpreter', 'latex', 'fontsize', 20)

figure
hold on
h2 = plot(n_dof, P_diff_fenics,'-k*', 'LineWidth', 2);
h1 = plot(xlim_vec, ones(1,length(xlim_vec))*p_diff_l,'--r', 'LineWidth', 2);
h1 = plot(xlim_vec, ones(1,length(xlim_vec))*p_diff_u,'--r', 'LineWidth', 2);
hold off
legend([h1,h2], {'Turek et al.', 'Fenics_IPCS'},'interpreter', ...
        'latex', 'fontsize', 16, 'Location','Best');
set(gca, 'XScale', 'log')
xlim(xlim_vec);
xlabel('ndofs', 'interpreter', 'latex', 'fontsize', 20)
ylabel('Pressure Difference', 'interpreter', 'latex', 'fontsize', 20)
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
