clc,clear

%% model parameters
alp = 0.5;
% sigma = 0.4; 
c = 1; % mobility
epsilon = 0.01; % interaction length

%% computational domain
T = 100;
a = 0;  b = 1;

%% mesh graded parameter
Nt = 5 * 10^(-1) * T;
%gam = 1;
gam = (2-alp)/alp; % 网格分级参数

%% time partition
T0 = min( [ 1/gam, T ]' );
T0 = T;
N0 = ceil( Nt/( T + 1 - 1/gam ) );
N0 = Nt;

% [0,T0]——graded mesh
t0 = T0 * ((0:1:N0)'/N0).^(gam);  % 时间分级网格
tau0 = diff( t0 );  % 分级网格步长

% [T0, T]——uniform mesh
N1 = Nt - N0;
%s_t = rand(N1,1);
s_t = ones(N1,1);
tau1 = (T - T0) * s_t/(sum(s_t));

tau = [tau0; tau1];
t = t0;
for k = 1:N1
    t(end+1,1) = t(end,1) + tau1( k, 1 );
end

%% spatial mesh partition
Nx = 128; 
hx = (b-a)/Nx;  % mesh size
x = (a+hx:hx:b)';  % mesh grid

%% stabilization parameter
S = 2; 

%% initial data
uu = load('u_intal_08.mat');
u_2 = uu.u_2;
u = reshape( u_2, Nx^2, 1 );
u_max(1,1) = max( abs( u(:,1) ) );

% %figure(1);
[X, Y] = meshgrid(x, x);
% % contourf(X, Y, u_2, 5)
% % colormap( 'jet' );
% % colorbar;
% 
% surf(X,Y,u_2)
% shading interp
% colormap( 'jet' );
% % colorbar
% % view([90, 90]);
% % axis off;
