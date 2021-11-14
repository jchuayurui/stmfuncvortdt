function [stmfunc] = solve_Poisson(vort,A,Nx,Ny)
%SOLVE_POISSON Summary of this function goes here
%   Detailed explanation goes here

% RHS of the Poisson equation. Note that the bc for streamfunction is zero.
b = -vort(:);

Ngs=(Nx-1)*(Ny-1);
u0=zeros(Ngs,1);
resobj=10^-4;

A=sparse(A);
stmfunc=iterativemethods(3,A,b,u0,resobj);

% When we solve the Poisson equation, we stack the 2D grid points into a 1D
% vector. But when we solve the time-dependent equation, we can work with a
% 2D mesh. So here, we reshape the 1D vector back to the 2D mesh
stmfunc = reshape(stmfunc,Nx-1,Ny-1);

end

