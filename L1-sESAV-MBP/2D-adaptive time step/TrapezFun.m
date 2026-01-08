function [ F ] = TrapezFun( hx, hy, u )
%% Numerical integration — Trapezoidal rule

%f(N+1,1) = f(1,1);

Nx = sqrt(length(u));

f = reshape(u, Nx, Nx);
ff = ( f(1:end-1,:) + f(2:end,:) )/2;
fff = ( ff(:,1:end-1) + ff(:,2:end) )/2;

F = hx * hy * sum( sum( fff, 1 ), 2 );

end

