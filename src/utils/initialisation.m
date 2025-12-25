function [u0] = initialisation(u,m,n)
%initialisation.m Sets the initial conditions of a Bessel vibration mode
%Inputs : 
%   u   : 3 dimensional matrix of the wave function
%   m   : Bessel azimuthal mode number (integer)
%   n   : Bessel radial mode number (integer)


global R;

[nr ntheta nt] = size(u);
dr = R/(nr-1);
dtheta = 2*pi/(ntheta-1);
r = 0:dr:R;
theta = 0:dtheta:2*pi;


%let's find lambda_(m,n)

f = @(x) besselj(m,x);
x0 = (n + m/2 - 1/4)*pi; 
jmn = fzero(f,x0);

lambda_mn = jmn / R;


u0 = u;
u0(:,:,1) = transpose(besselj(m, lambda_mn*r))*cos(m*theta);

 


end