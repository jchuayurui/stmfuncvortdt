function [vortnew] = advance_vort(stmfunc,vort,Nx,Ny,dx,dy,dt,Re,t)
%SOLVE_ Summary of this function goes here
%   Detailed explanation goes here

RHS = assembleRHS(Nx,Ny,stmfunc,vort,Re,dx,dy,t);

vortnew = vort + dt*RHS;
                                
end









