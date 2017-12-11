% Felix Newberry

% Compare FSI fluid domain to Lid driven cavity

clear all
close all
clc



% Ghiea et al
% Results for u velocity along Vertical Line through Geometric Center of
% Cavityu_out
u_ghia_100 = [1, 0.84123, 0.78871, 0.73722, 0.68717, 0.23151, 0.00332, ...
    -0.13641, -0.20581, -0.21090, -0.15662, -0.10150, -0.06434, ...
    -0.04775, -0.04192, -0.03717, 0];
y_ghia =[1, 0.9766, 0.9688, 0.9609, 0.9531, 0.8516, 0.7344, 0.6172, ...
    0.5, 0.4531, 0.2813, 0.1719, 0.1016, 0.0703, 0.0625, 0.0547, 0];

v_ghia_100 = [0, -0.05906, -0.07391, -0.08864, -0.10313, -0.16914, ...
    -0.22445, -0.24533, 0.05454, 0.17527, 0.17507, 0.16077, 0.12317, ...
    0.10890, 0.10091, 0.09233, 0];
x_ghia = [1, 0.9688, 0.9609, 0.9531, 0.9453, 0.9063, 0.8594, 0.8047, ...
    0.5, 0.2344, 0.2266, 0.1563, 0.0938, 0.0781, 0.0703, 0.0625, 0];


% FSI results
% velocities at y = 0.5 in y direction as x changes
load('u_u.mat')
load('u_v.mat')

% lid driven cavity chorin results
load('u_u_chorin.mat')
load('u_v_chorin.mat')

x = u_v(:,1); 
v_sim = u_v(:,2); 

y = u_u(:,1)-0.2;
u_sim = u_u(:,2); 

x_c = u_v_chorin(:,1); 
v_sim_c = u_v_chorin(:,2); 

y_c = u_u_chorin(:,1);
u_sim_c = u_u_chorin(:,2); 

figure
subplot(2,1,1)
hold on
p1 = plot(u_sim, y, 'c-', 'LineWidth', 2);
p2 = plot(u_sim_c, y_c, 'k-', 'LineWidth', 1);
p3 = plot(u_ghia_100, y_ghia,'ro', 'LineWidth', 2);
xlabel('u velocity', 'interpreter', 'latex', 'fontsize', 20)
ylabel('y', 'interpreter', 'latex', 'fontsize', 20)
grid on
hold off

subplot(2,1,2)
hold on
plot(x,v_sim, 'c-', 'LineWidth', 2)
plot(x_c,v_sim_c, 'k-', 'LineWidth', 1)
plot(x_ghia, v_ghia_100,'ro', 'LineWidth', 2)
xlabel('x', 'interpreter', 'latex', 'fontsize', 20)
ylabel('v velocity', 'interpreter', 'latex', 'fontsize', 20)
legend([p1, p2, p3],{'fenics','chorin ldc','ghia et al'},'interpreter', 'latex', 'fontsize', 16)

grid on