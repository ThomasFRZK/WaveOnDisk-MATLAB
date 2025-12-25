
%-------------------------------------------%
%           DISCRETISATION of EPD           %
%-------------------------------------------%

%{
    In this part we will focus on the discretisation of the D'Alembert EPD.
Let u be the wave fonction solution of our EPD, and let us call dr, dtheta
and dt the increment steps for variables r, theta and t respectively. The
radius of the disk and the time of simulation are written as R and T.

%}


%--------------------%
%   INITIALISATION   %
%--------------------%

global R T c;
R = 1;
T = 2;
c = 1;

%The increment steps : 

dr = 0.05;
dtheta = 2*pi/100;
dt = 0.5*dr*dtheta;

% As u(r,theta,t) is a function depending of three variables, let's build a
% 3D matrix U.

%Let us define the size of each dimension : 

Nt = round(T/dt);
Nr = round(R/dr);
Ntheta = round(2*pi/dtheta);


U = zeros(Nr+1,Ntheta+1,Nt); %Add 1 to Nr to take into account r = 0.

%First, let's highlight different disk normal modes of vibration : 

%Here, one could choose the Bessel vibrationnal modes of the disk

M_bessel = 0; 
N_bessel = 2;


U = initialisation(U, M_bessel, N_bessel);

%Then, we calculate the Taylor expansion to second order at t = delta t

U(:,:,2) = U(:,:,1)+dt^2/2 * c^2 * laplacian(U(:,:,1));

for n = 3:Nt
    U(:,:,n) = c^2*dt^2*laplacian(U(:,:,n-1))+2*U(:,:,(n-1))-U(:,:,(n-2));
end


%-------------------------------------%
%       VALIDATION OF THE MODEL       %
%-------------------------------------%

%{
    Now let us confirm that the D'Alembert equation leads correctly to the
actual Bessel vibration modes. First, we will implement the mean square
error (MSE) and then we will check under what conditions our model may
diverge.

%}

%Define the matrix of MSEs : 
Mses = zeros(1,Nt);

%Now let's define the real solution matrix

U_real = BesselSolution(M_bessel,N_bessel, Nr, Ntheta, Nt);

for n = 1:Nt
    Mses(n) = MSE(U(:,:,n), U_real(:,:,n));
end



%-----------------------------------------%
%       SUPPERPOSITION OF IMPULSES        %
%-----------------------------------------%

%{
    Let us now play with different linear combinations of Bessel functions
 (which is still a solution of the PDE) to see how behaves the disk.

%}

%Let's define the amplitude coefficients we want to use
C = [0.1 0.4 0.25 0.25];

%Now, we configure the differents Bessel modes of our waves

M_modes = [ 0 1 2 3];
N_modes = [ 2 3 2 0];

%Then, the solution matrix is 
V = pulsesCombination(C, M_modes,N_modes, Nr, Ntheta, Nt);


%---------------------------%
%       PLOT SECTION        %
%---------------------------%


%Let's convert our cylindrical coordonates matrix into Cartesian
%coordinates.

r = (0:Nr)'*dr;
theta = (0:Ntheta)*dtheta;

% Cylindrical grid
[R_grid, Theta_grid] = meshgrid(r, theta);
R_grid = R_grid';
Theta_grid = Theta_grid';

% Conversion cylindrical -> cartesian
X = R_grid .* cos(Theta_grid);
Y = R_grid .* sin(Theta_grid);

% --- Animated vibration mode --- %

figure(1)
hSurf = surf(X, Y, U(:,:,1));
axis([-1 1 -1 1 -1 1]);
shading interp
colormap jet
colorbar
xlabel('X')
ylabel('Y')
zlabel('u(r,\theta)')
title('Animated wave on the disk')

%Let us animate the figure over time
for t = 2:size(U,3)
    hSurf.ZData = U(:,:,t);  
    drawnow;
    pause(0.05);
end

% --- MSE plot --- %

figure(2);
plot((0:Nt-1)*dt, Mses,'o');
xlabel("Time");
ylabel("MSE");
title("Mean Square Error in time of the whole wave")


% --- Plot of the frequency gap --- %

figure(3);
t = 0:dt:(T-dt);
r0 = 10;
theta0 = 14;

u_num0 = reshape(U(r0,theta0,:), size(t));
u_real0 = reshape(U_real(r0,theta0,:), size(t));

plot(t,u_num0, t, u_real0);
legend("Numerical wave", "Real wave")
title(sprintf("disk altitude at r = %d and theta = %d in time", r0*dr, theta0*dtheta));



ynum = fft(u_num0);
yreal = fft(u_real0);

figure(4)
plot(t, abs(ynum), t, abs(yreal));
legend("Numerical wave", "Real wave")

title("Fast Fourier Transform of real and numerical waves")

% --- Linear Combination of modes animation -- %

figure(5)
hSurf = surf(X, Y, V(:,:,1));
axis([-1 1 -1 1 -1 1]);
shading interp
colormap jet
colorbar
xlabel('X')
ylabel('Y')
zlabel('u(r,\theta)')
title('Animated wave combination on the disk')


for t = 2:size(V,3)
    hSurf.ZData = V(:,:,t);  
    drawnow;
    pause(0.05);
end
