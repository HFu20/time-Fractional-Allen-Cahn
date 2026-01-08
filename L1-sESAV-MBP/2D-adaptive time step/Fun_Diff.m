function [ FF, F, dF ] = Fun_Diff( u )
% information related to the nonlinear term

FF = 1/4 * ( u.^2 - 1 ).^2;
F = u - u.^3;
dF = 1 - 3 * u.^2;

end

