function [RHS] = assembleRHS(Nx,Ny,stmfunc,vort,Re,dx,dy,t)
%ASSEMBLERHS Summary of this function goes here
%   Detailed explanation goes here

% the boundary conditions at the boundaries
U_south = zeros(Nx-1,1);
U_north = sin(pi*t*10)*ones(Nx-1,1);
U_west  = zeros(Ny-1,1);
U_east  = zeros(Ny-1,1);

nu=1/Re;
dx2=dx*dx;
dy2=dy*dy;
RHS=zeros(Nx-1,Ny-1);

% interior
for j=2:Ny-2
    for i=2:Nx-2
        fac1 = -(stmfunc(i,j+1) - stmfunc(i,j-1))/2/dy;
        fac2 =  (stmfunc(i+1,j) - stmfunc(i-1,j))/2/dx;

        RHS(i,j) =  fac1*( vort(i+1,j)-vort(i-1,j) )/2/dx + ...
                    fac2*( vort(i,j+1)-vort(i,j-1) )/2/dy + ...
                    nu*( ( vort(i+1,j)-2*vort(i,j)+vort(i-1,j) )/dx2 + ... 
                         ( vort(i,j+1)-2*vort(i,j)+vort(i,j-1) )/dy2 );                     
    end
end

% south
j=1;

for i=2:Nx-2
    vortsouthbc = ( 0 - stmfunc(i,j) + U_south(i)*dy)/(0.5*dy2);
        
    fac1 = -(stmfunc(i,j+1) - 0             )/2/dy;
    fac2 =  (stmfunc(i+1,j) - stmfunc(i-1,j))/2/dx;

    RHS(i,j) =  fac1*( vort(i+1,j)-vort(i-1,j) )/2/dx + ...
                fac2*( vort(i,j+1)-vortsouthbc )/2/dy + ...
                nu*( ( vort(i+1,j)-2*vort(i,j)+vort(i-1,j) )/dx2 + ... 
                     ( vort(i,j+1)-2*vort(i,j)+vortsouthbc )/dy2 ); 
end

% North
j=Ny-1;
for i=2:Nx-2
    vortnorthbc = ( 0 - stmfunc(i,j) - U_north(i)*dy)/(0.5*dy2);
        
    fac1 = -(0              - stmfunc(i,j-1))/2/dy;
    fac2 =  (stmfunc(i+1,j) - stmfunc(i-1,j))/2/dx;

    RHS(i,j) =  fac1*( vort(i+1,j)-vort(i-1,j) )/2/dx + ...
                fac2*( vortnorthbc-vort(i,j-1) )/2/dy + ...
                nu*( ( vort(i+1,j)-2*vort(i,j)+vort(i-1,j) )/dx2 + ... 
                     ( vortnorthbc-2*vort(i,j)+vort(i,j-1) )/dy2); 
end

% West
i=1;
for j=2:Ny-2
    vortwestbc = ( 0 - stmfunc(i,j) - U_west(j)*dx)/(0.5*dx2);
    
    fac1 = -(stmfunc(i,j+1) - stmfunc(i,j-1))/2/dy;
    fac2 =  (stmfunc(i+1,j) - 0             )/2/dx;
 
    RHS(i,j) =  fac1*( vort(i+1,j)-vortwestbc  )/2/dx + ...
                fac2*( vort(i,j+1)-vort(i,j-1) )/2/dy + ...
                nu*( ( vort(i+1,j)-2*vort(i,j)+vortwestbc  )/dx2 + ... 
                     ( vort(i,j+1)-2*vort(i,j)+vort(i,j-1) )/dy2 ); 
end

% East
i=Nx-1;
for j=2:Ny-2
    vorteastbc = ( 0 - stmfunc(i,j) + U_east(j)*dx)/(0.5*dx2);
    
    fac1 = -(stmfunc(i,j+1) - stmfunc(i,j-1))/2/dy;
    fac2 =  (0              - stmfunc(i-1,j))/2/dx;
 
    RHS(i,j) =  fac1*( vorteastbc -vort(i-1,j) )/2/dx + ...
                fac2*( vort(i,j+1)-vort(i,j-1) )/2/dy + ...
                nu*( ( vorteastbc -2*vort(i,j)+vort(i-1,j) )/dx2 + ... 
                     ( vort(i,j+1)-2*vort(i,j)+vort(i,j-1) )/dy2 ); 
end

% South-west
i=1;j=1;

    vortsouthbc = ( 0 - stmfunc(i,j) + U_south(i)*dy)/(0.5*dy2);
    vortwestbc  = ( 0 - stmfunc(i,j) + U_west(j) *dx)/(0.5*dx2);
    
    fac1 = -(stmfunc(i,j+1) - 0 )/2/dy;
    fac2 =  (stmfunc(i+1,j) - 0 )/2/dx;

    RHS(i,j) =  fac1*( vort(i+1,j)-vortwestbc  )/2/dx + ...
                fac2*( vort(i,j+1)-vortsouthbc )/2/dy + ...
                nu*( ( vort(i+1,j)-2*vort(i,j)+vortwestbc  )/dx2 + ... 
                     ( vort(i,j+1)-2*vort(i,j)+vortsouthbc )/dy2 ); 
                    
% South-east          
i=Nx-1;j=1;

    vortsouthbc = ( 0 - stmfunc(i,j) + U_south(i)*dy)/(0.5*dy2);
    vorteastbc  = ( 0 - stmfunc(i,j) - U_east(j) *dx)/(0.5*dx2);
    
    fac1 = -(stmfunc(i,j+1) - 0)/2/dy;
    fac2 =  (0 - stmfunc(i-1,j))/2/dx;

    RHS(i,j) =  fac1*( vorteastbc -vort(i-1,j) )/2/dx + ...
                fac2*( vort(i,j+1)-vortsouthbc )/2/dy + ...
                nu*( ( vorteastbc -2*vort(i,j)+vort(i-1,j) )/dx2 + ... 
                     ( vort(i,j+1)-2*vort(i,j)+vortsouthbc )/dy2 ); 
                    
% North-east         
i=Nx-1;j=Ny-1;

    vortnorthbc = ( 0 - stmfunc(i,j) - U_north(i)*dy)/(0.5*dy2);
    vorteastbc  = ( 0 - stmfunc(i,j) - U_east(j) *dx)/(0.5*dx2);
        
    fac1 = -(0 - stmfunc(i,j-1))/2/dy;
    fac2 =  (0 - stmfunc(i-1,j))/2/dx;

    RHS(i,j) =  fac1*( vorteastbc -vort(i-1,j) )/2/dx + ...
                fac2*( vortnorthbc-vort(i,j-1) )/2/dy + ...
                nu*( ( vorteastbc -2*vort(i,j)+vort(i-1,j) )/dx2 + ... 
                     ( vortnorthbc-2*vort(i,j)+vort(i,j-1) )/dy2 );

% North-west          
i=1;j=Ny-1;

    vortnorthbc = ( 0 - stmfunc(i,j) - U_north(i)*dy)/(0.5*dy2);
    vortwestbc  = ( 0 - stmfunc(i,j) + U_west(j) *dx)/(0.5*dx2);
        
    fac1 = -(0 - stmfunc(i,j-1))/2/dy;
    fac2 =  (stmfunc(i+1,j) - 0)/2/dx;

    RHS(i,j) =  fac1*( vort(i+1,j)-vortwestbc  )/2/dx + ...
                fac2*( vortnorthbc-vort(i,j-1) )/2/dy + ...
                nu*( ( vort(i+1,j)-2*vort(i,j)+vortwestbc  )/dx2 + ... 
                     ( vortnorthbc-2*vort(i,j)+vort(i,j-1) )/dy2 ); 
                                                       
end

