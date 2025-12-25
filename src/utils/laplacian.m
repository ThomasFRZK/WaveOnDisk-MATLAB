function [Lap_u] = laplacian(u)
%laplacian.m calculates the laplacien in spherical coordonates for the
%matrix u.
%Input :
%   u   : 2D matrix of wave amplitude values.

global R;

[nr ntheta] = size(u);
dr = R/(nr-1);
dtheta = 2*pi/(ntheta-1);

Lap_u = zeros(size(u));

Lap_u(1,:) = 2*(mean(u(2,:))-u(1,:))/dr^2;

for n = 2:nr-1
    r = n*dr;
    for m = 1:ntheta

        dr_u = (u(n+1,m)-u(n-1,m))/(2*dr);
        drr_u = (u(n+1,m)-2*u(n,m)+u(n-1,m))/dr^2;
        dthetaa_u = (u(n,mod(m,ntheta)+1) - 2*u(n,m) + u(n,mod(m-2,ntheta)+1))/dtheta^2;
        
        Lap_u(n,m) = drr_u + 1/r *dr_u + 1/r^2 * dthetaa_u;
    
    end
end


end