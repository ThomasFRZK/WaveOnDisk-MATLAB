function [u] = pulsesCombination(coeffs, m_array, n_array, nr, ntheta, nt)
%pulsesCombination.m gives the solution matrix of a combination of Bessel
%functions on the disk.
%INPUTS : 
%   coeffs  : Array of amplitude of each Bessel wave
%   m_array : Array of Bessel azimuthal numbers
%   n_array : Array of Bessel radial numbers
%   nr      : Number of radial knots
%   ntheta  : Number of angular knots
%   nt      : Number of time steps

u = 0;
N = length(coeffs);
for i=1:N
    u = coeffs(i)*BesselSolution(m_array(i), n_array(i), nr, ntheta, nt)+u;
end 

end