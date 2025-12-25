function [u] = BesselSolution(m,n,nr,ntheta,nt)
%BesselSolution.m creates a nr x ntheta x nt array of the real D'Alembert
%equation solution for vibrational modes (m,n)

global R T c;



r = linspace(0,R,nr+1);
theta = linspace(0,2*pi,ntheta+1);
t = linspace(0,T,nt);

%let's find lambda_(m,n) 
f  = @(x) besselj(m,x);
x0 = (n + m/2 - 1/4)*pi;
jmn = fzero(f,x0);

lambda_mn = jmn / R;
omega_mn  = c * lambda_mn;

Rpart = repmat(transpose(besselj(m, lambda_mn * r)), [1 (ntheta+1) nt]);        
Thetapart = repmat(cos(m * theta), [(nr+1) 1 nt]) ;                          
Tpart = repmat(reshape(cos(omega_mn * t), 1, 1, nt), (nr+1), (ntheta+1),1);
u = Rpart .* Thetapart .* Tpart;

end
