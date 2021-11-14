function [vortnew] = advance_vortimp(stmfunc,vort,Nx,Ny,dx,dy,dt,Re,t)
%SOLVE_ Summary of this function goes here
%   Detailed explanation goes here

[Aimp,bimp]=assembleAimp(stmfunc,Nx,Ny,dx,dy,Re,dt,t);
vort = reshape(vort,(Nx-1)*(Ny-1),1);
vect = (vort/dt) + bimp;
vortnew = inv(Aimp)*vect;
                                
end









