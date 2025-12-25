function [err_] = MSE(u_num, u_real)
%MSE.m calculates the mean square error between two arrays
%INPUTS : (with the same dimensions)
%   u_num   : array corresponding to the numerical solution
%   u_real  : array corresponding to the real solution

    [nr ntheta] = size(u_num);

    err_ = mean((u_num - u_real).^2, "all");

end