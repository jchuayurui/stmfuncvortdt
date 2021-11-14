function [A,b]=assembleAimp(stmfunc,Nx,Ny,dx,dy,Re,dt,t)
b=zeros((Nx-1)*(Ny-1),1);
dx2=dx*dx;
dy2=dy*dy;
U_south = zeros(Nx-1,1);
U_north = sin(pi*t*10)*ones(Nx-1,1);
U_west  = zeros(Ny-1,1);
U_east  = zeros(Ny-1,1);
A=zeros((Nx-1)*(Ny-1),(Nx-1)*(Ny-1));

% interior
for j=2:Ny-2
    for i=2:Nx-2
        fac1 = -(stmfunc(i,j+1) - stmfunc(i,j-1))/2/dy;
        fac2 =  (stmfunc(i+1,j) - stmfunc(i-1,j))/2/dx;
        po=i+(j-1)*(Nx-1);
        A(po,po)= (1/dt)+((2*dy2+2*dx2)/(Re*dx2*dy2)); 
        A(po,po+1)=(-Re*fac1*dx-2)/(2*Re*dx2); 
        A(po,po-1)=(Re*fac1*dx-2)/(2*Re*dx2);
        A(po,po-(Nx-1))=(Re*fac2*dy-2)/(2*Re*dy2);
        A(po,po+(Nx-1))=(-Re*fac2*dy-2)/(2*Re*dy2);
    end
end

% South
j=1;
for i=2:Nx-2
	vortsouthbc = ( 0 - stmfunc(i,j) + U_south(i)*dy)/(0.5*dy2);
        
    fac1 = -(stmfunc(i,j+1) - 0             )/2/dy;
    fac2 =  (stmfunc(i+1,j) - stmfunc(i-1,j))/2/dx;
        po=i+(j-1)*(Nx-1);
        A(po,po)= (1/dt)+((2*dy2+2*dx2)/(Re*dx2*dy2)); 
        A(po,po+1)=(-Re*fac1*dx-2)/(2*Re*dx2); 
        A(po,po-1)=(Re*fac1*dx-2)/(2*Re*dx2);
      %  A(po,po-(Nx-1))=(Re*fac2*dy-2)/(2*Re*dy2);
        A(po,po+(Nx-1))=(-Re*fac2*dy-2)/(2*Re*dy2);
        b(po) = -((Re*fac2*dy-2)/(2*Re*dy2))* vortsouthbc;
end

% North
j=Ny-1;
for i=2:Nx-2
	vortnorthbc = ( 0 - stmfunc(i,j) - U_north(i)*dy)/(0.5*dy2);
        
    fac1 = -(0              - stmfunc(i,j-1))/2/dy;
    fac2 =  (stmfunc(i+1,j) - stmfunc(i-1,j))/2/dx;
        po=i+(j-1)*(Nx-1);
        A(po,po)= (1/dt)+((2*dy2+2*dx2)/(Re*dx2*dy2)); 
        A(po,po+1)=(-Re*fac1*dx-2)/(2*Re*dx2); 
        A(po,po-1)=(Re*fac1*dx-2)/(2*Re*dx2);
        A(po,po-(Nx-1))=(Re*fac2*dy-2)/(2*Re*dy2);
        %A(po,po+(Nx-1))=(-Re*fac2*dy-2)/(2*Re*dy2);
        b(po) = -((-Re*fac2*dy-2)/(2*Re*dy2))* vortnorthbc;
end

% West
i=1;
for j=2:Ny-2
	vortwestbc = ( 0 - stmfunc(i,j) - U_west(j)*dx)/(0.5*dx2);
    
    fac1 = -(stmfunc(i,j+1) - stmfunc(i,j-1))/2/dy;
    fac2 =  (stmfunc(i+1,j) - 0             )/2/dx;
        po=i+(j-1)*(Nx-1);
        A(po,po)= (1/dt)+((2*dy2+2*dx2)/(Re*dx2*dy2)); 
        A(po,po+1)=(-Re*fac1*dx-2)/(2*Re*dx2); 
       % A(po,po-1)=(Re*fac1*dx-2)/(2*Re*dx2);
        A(po,po-(Nx-1))=(Re*fac2*dy-2)/(2*Re*dy2);
        A(po,po+(Nx-1))=(-Re*fac2*dy-2)/(2*Re*dy2);
        b(po) = -((Re*fac1*dx-2)/(2*Re*dx2))* vortwestbc;
end

% East
i=Nx-1;
for j=2:Ny-2
	vorteastbc = ( 0 - stmfunc(i,j) + U_east(j)*dx)/(0.5*dx2);
    
    fac1 = -(stmfunc(i,j+1) - stmfunc(i,j-1))/2/dy;
    fac2 =  (0              - stmfunc(i-1,j))/2/dx;
        po=i+(j-1)*(Nx-1);
        A(po,po)= (1/dt)+((2*dy2+2*dx2)/(Re*dx2*dy2));  
       % A(po,po+1)=(-Re*fac1*dx-2)/(2*Re*dx2); 
        A(po,po-1)=(Re*fac1*dx-2)/(2*Re*dx2);
        A(po,po-(Nx-1))=(Re*fac2*dy-2)/(2*Re*dy2);
        A(po,po+(Nx-1))=(-Re*fac2*dy-2)/(2*Re*dy2);
        b(po) = -((-Re*fac1*dx-2)/(2*Re*dx2))* vorteastbc;
end

% South-west
i=1;j=1;
	vortsouthbc = ( 0 - stmfunc(i,j) + U_south(i)*dy)/(0.5*dy2);
    vortwestbc  = ( 0 - stmfunc(i,j) + U_west(j) *dx)/(0.5*dx2);
    
    fac1 = -(stmfunc(i,j+1) - 0 )/2/dy;
    fac2 =  (stmfunc(i+1,j) - 0 )/2/dx;
        po=i+(j-1)*(Nx-1);
        A(po,po)= (1/dt)+((2*dy2+2*dx2)/(Re*dx2*dy2)); 
        A(po,po+1)=(-Re*fac1*dx-2)/(2*Re*dx2); 
       % A(po,po-1)=(Re*fac1*dx-2)/(2*Re*dx2);
       % A(po,po-(Nx-1))=(Re*fac2*dy-2)/(2*Re*dy2);
        A(po,po+(Nx-1))=(-Re*fac2*dy-2)/(2*Re*dy2);
        b(po) = (-((Re*fac1*dx-2)/(2*Re*dx2))* vortwestbc)-(((Re*fac2*dy-2)/(2*Re*dy2))* vortsouthbc);
% South-east          
i=Nx-1;j=1;
	vortsouthbc = ( 0 - stmfunc(i,j) + U_south(i)*dy)/(0.5*dy2);
    vorteastbc  = ( 0 - stmfunc(i,j) - U_east(j) *dx)/(0.5*dx2);
    
    fac1 = -(stmfunc(i,j+1) - 0)/2/dy;
    fac2 =  (0 - stmfunc(i-1,j))/2/dx;
        po=i+(j-1)*(Nx-1);
        A(po,po)= (1/dt)+((2*dy2+2*dx2)/(Re*dx2*dy2)); 
     %   A(po,po+1)=(-Re*fac1*dx-2)/(2*Re*dx2); 
        A(po,po-1)=(Re*fac1*dx-2)/(2*Re*dx2);
       % A(po,po-(Nx-1))=(Re*fac2*dy-2)/(2*Re*dy2);
        A(po,po+(Nx-1))=(-Re*fac2*dy-2)/(2*Re*dy2);
		b(po) = (-((Re*fac2*dy-2)/(2*Re*dy2))* vortsouthbc) - (((-Re*fac1*dx-2)/(2*Re*dx2))* vorteastbc);
% North-east         
i=Nx-1;j=Ny-1;
	vortnorthbc = ( 0 - stmfunc(i,j) - U_north(i)*dy)/(0.5*dy2);
    vorteastbc  = ( 0 - stmfunc(i,j) - U_east(j) *dx)/(0.5*dx2);
        
    fac1 = -(0 - stmfunc(i,j-1))/2/dy;
    fac2 =  (0 - stmfunc(i-1,j))/2/dx;
        po=i+(j-1)*(Nx-1);
        A(po,po)= (1/dt)+((2*dy2+2*dx2)/(Re*dx2*dy2)); 
     %   A(po,po+1)=(-Re*fac1*dx-2)/(2*Re*dx2); 
        A(po,po-1)=(Re*fac1*dx-2)/(2*Re*dx2);
        A(po,po-(Nx-1))=(Re*fac2*dy-2)/(2*Re*dy2);
       % A(po,po+(Nx-1))=(-Re*fac2*dy-2)/(2*Re*dy2);
       b(po) = (-((-Re*fac2*dy-2)/(2*Re*dy2))* vortnorthbc)- (((-Re*fac1*dx-2)/(2*Re*dx2))* vorteastbc);
% North-west          
i=1;j=Ny-1;
	vortnorthbc = ( 0 - stmfunc(i,j) - U_north(i)*dy)/(0.5*dy2);
    vortwestbc  = ( 0 - stmfunc(i,j) + U_west(j) *dx)/(0.5*dx2);
        
    fac1 = -(0 - stmfunc(i,j-1))/2/dy;
    fac2 =  (stmfunc(i+1,j) - 0)/2/dx;
        po=i+(j-1)*(Nx-1);
        A(po,po)= (1/dt)+((2*dy2+2*dx2)/(Re*dx2*dy2)); 
        A(po,po+1)=(-Re*fac1*dx-2)/(2*Re*dx2); 
     %   A(po,po-1)=(Re*fac1*dx-2)/(2*Re*dx2);
        A(po,po-(Nx-1))=(Re*fac2*dy-2)/(2*Re*dy2);
       % A(po,po+(Nx-1))=(-Re*fac2*dy-2)/(2*Re*dy2);
        b(po) = (-((-Re*fac2*dy-2)/(2*Re*dy2))* vortnorthbc)-(((Re*fac1*dx-2)/(2*Re*dx2))* vortwestbc);

end
