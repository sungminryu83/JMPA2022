clc
clear all

% The y coordinate of the test interval
v_coor = pi/6; 

% The number of sample points
n_point = 500; 

% The x coordiantes of the start and end points of the test interval
u_i = 0;  u_e = pi;   

% Generating a vector for the test interval
u_vec_tmp = linspace(u_i, u_e, n_point); 
u_vec = [u_vec_tmp(2:end-1)];

% A parameterization for the sphere with radius 1
syms u v real;
f = sin(u)*cos(v);
g = sin(u)*sin(v);
h = cos(u);
X = [f g h];

% The first and second order derivatives of r
Xu = diff(X, u, 1);
Xv = diff(X, v, 1);
Xuu = diff(X, u, 2);
Xuv = diff(X, u, v);
Xvv = diff(X, v, 2);

% The coefficients of the first fundamental form
E = dot(Xu, Xu);
F = dot(Xu, Xv);
G = dot(Xv, Xv);

% The surface normal vector and the normalization of this vector
Nor = cross(Xu,Xv);
n = Nor/norm(Nor,2);

% The coefficients of the second fundamental form
L = dot(Xuu, n);
M = dot(Xuv, n);
N = dot(Xvv, n);

% Computing the Gaussian curvature
K = (L*N-M^2)/(E*G-F^2); 

% Computing the geodesic curvatures over the coordinate curves
kg1 = dot(Xuu,cross(n,Xu))/(norm(Xu,2))^3;
kg2 = dot(Xvv,cross(n,Xv))/(norm(Xv,2))^3;

% Computing Fa and Fb
Fa = kg1*sqrt(E);
Fb = kg2*sqrt(G);

% Computing the anlge of intersection 
phi = acos(F/sqrt(E*G));

% Computing I, II, III, and I+II+III
I = simplify(K*norm(Nor,2));
II = simplify(diff(Fb,u)-diff(Fa,v));
III = simplify(-1*diff(phi,u,v));
net = I+II+III;
net = simplify(net);

% Setting functions for I, II, III, and I+II+III
syms I_f(u,v) II_f(u,v) III_f(u,v) net_f(u,v);
I_f(u,v) = I;
II_f(u,v) = II; 
III_f(u,v) = III; 
net_f(u,v) = net;

% Computing the values of the functions of I, II, III, and I+II+III for the test interval
I_vec = double(I_f(u_vec,v_coor));
II_vec = double(II_f(u_vec,v_coor));
III_vec = double(III_f(u_vec,v_coor));
net_vec = I_vec + II_vec + III_vec;

% Setting parameters for figure configuration
xlabelfontsize_value = 18;
ylabelfontsize_value = 18;
ticklabelfontsize_value = 18;
legendfontsize_value = 13;
linewidth_value = 1.3;
xtick_interval = 1.0;
ytick_interval = 1.5;

% Plotting 
figure(1)
plot(u_vec, I_vec, '-- k', 'LineWidth',linewidth_value); hold on;
plot(u_vec, II_vec, ': k', 'LineWidth',linewidth_value);
plot(u_vec, III_vec, '- k', 'LineWidth',linewidth_value);
%plot(u_vec, net_vec, '- k', 'LineWidth',linewidth_value);

x_i = 0; x_e = pi;
y_i = -1.5; y_e = 2.5;

xlim([x_i x_e]); 
ylim([y_i y_e]);

set(gca,'Xtick',x_i:xtick_interval:x_e,'FontSize', ticklabelfontsize_value,'FontName','Times');
set(gca,'Ytick',y_i:ytick_interval:y_e,'FontSize', ticklabelfontsize_value,'FontName','Times');

xlabel('$u$','Interpreter','Latex','FontSize', xlabelfontsize_value); 
ylabel('Budgets of (2.7)','Interpreter','Latex','FontSize', ylabelfontsize_value);
legend('$I_{K}$','$I_{\kappa_{g}}$','$I_{\phi}$','Interpreter','Latex', 'FontSize',legendfontsize_value);
