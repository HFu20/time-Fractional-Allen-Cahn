clc,clear

%% model parameters
alp = 0.5;
% sigma = 0.4; 
c = 1; % mobility
epsilon = 0.01; % interaction length

%% computational domain
T = 5;
a = 0;  b = 1;

%% mesh graded parameter
% Nt = 5 * 10^(-1) * T;
% gam = 1;
gam = (2-alp)/alp; % 网格分级参数

%% time partition
T0 = 0.5;
N0 = 30;

% [0,T0]——graded mesh
t = T0 * ((0:1:N0)'/N0).^(gam);  % 时间分级网格
tau = diff( t );  % 分级网格步长

% [T0,T]——parameters in adaptive time-stepping method
tau_min = 0.02;
tau_max = 2;
bet = 10^7;
Nt = 10^13;

%% spatial mesh partition
Nx = 128;  % 空间网格剖分次数
hx = (b-a)/Nx;  % mesh size
x = (a+hx:hx:b)';  % mesh grid

%% stabilization parameter
S = 2; 

%% initial data
uu = load('u_intal_08.mat');
u_2 = uu.u_2;
u = reshape( u_2, Nx^2, 1 );
u_max(1,1) = max( abs( u(:,1) ) );

%figure(1);
[X, Y] = meshgrid(x, x);
% contourf(X, Y, u_2, 5)
% colormap( 'jet' );
% colorbar;

% surf(X,Y,u_2)
% shading interp
% colormap( 'jet' );
% colorbar
% view([90, 90]);
% axis off;
