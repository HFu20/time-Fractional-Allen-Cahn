clc,clear

%% 模型参数
alp = 0.9;
c = 1; % 迁移率
epsilon = 0.01; % 界面宽度

%% 离散点
sig = 1 - alp/2;

%% 求解区域
T = 100;
a = 0;  b = 1;

%% 时间网格剖分参数
% gam = 1;
gam = 2/alp; % 网格分级参数

%% 时间网格剖分
T0 = 1/3;
N0 = 1;
% T0 = 0;
% N0 = 0;

% [0,T0]——分级网格
t = T0 * ((0:1:N0)'/N0).^(gam);  % 时间分级网格
tau = diff( t );  % 分级网格步长

% [T0,T]——自适应网格
tau_min = 1/30;
tau_max = 1/3;
bet = 10^7;
Nt = 10^13;

for k = 1:N0-1
    rho(k,1) = tau(k,1)/tau(k+1,1);
end

%% 空间网格剖分
Nx = 128;  % 空间网格剖分次数
hx = (b-a)/Nx;  % mesh size
x = (a+hx:hx:b)';  % mesh grid

%% 稳定化
S = 8.02; % 稳定化常数

%% 问题初值
uu = load('u_intal_08.mat');
u_2 = uu.u_2;
u = reshape( u_2, Nx^2, 1 );
u_max(1,1) = max( abs( u(:,1) ) );

[X, Y] = meshgrid(x, x);
% figure(1);
% surf(X,Y,u_2)
% shading interp
% colormap( 'jet' );
% colorbar
% view([90, 90]);
% axis off;