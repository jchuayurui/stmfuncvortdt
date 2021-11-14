function [u,v] = get_uv(stmfunc,Nx,Ny,dx,dy,t)
%GET_UV Summary of this function goes here
%   Detailed explanation goes here

stmfunc=[0             zeros(1,Ny-1)   0
         zeros(Nx-1,1) stmfunc         zeros(Nx-1,1)
         0             zeros(1,Ny-1)   0 ]; 
     
u =  (stmfunc(2:end-1,3:end) - stmfunc(2:end-1,1:end-2))/2/dy;
v = -(stmfunc(3:end,2:end-1) - stmfunc(1:end-2,2:end-1))/2/dx;

u=[0             zeros(1,Ny-1)  0
   zeros(Nx-1,1) u              (sin(pi*t*10)*ones(Nx-1,1))
   0             zeros(1,Ny-1)  0 ]; 

v=[0             zeros(1,Ny-1)  0
   zeros(Nx-1,1) v              zeros(Nx-1,1)
   0             zeros(1,Ny-1)  0 ]; 

end

