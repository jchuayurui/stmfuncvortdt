function [A]=assembleA(Nx,Ny,dx,dy)

dx2=dx*dx;
dy2=dy*dy;

A=zeros((Nx-1)*(Ny-1),(Nx-1)*(Ny-1));

% interior
for j=2:Ny-2
    for i=2:Nx-2
        po=i+(j-1)*(Nx-1);
        A(po,po)= (-2/dx2)-(2/dy2); 
        A(po,po+1)=1/dx2; 
        A(po,po-1)=1/dx2;
        A(po,po-(Nx-1))=1/dy2;
        A(po,po+(Nx-1))=1/dy2;
    end
end

% South
j=1;
for i=2:Nx-2
        po=i+(j-1)*(Nx-1);
        A(po,po)= (-2/dx2)-(2/dy2); 
        A(po,po+1)=1/dx2; 
        A(po,po-1)=1/dx2;
        %A(po,po-(Nx-1))=1/dy2;
        A(po,po+(Nx-1))=1/dy2;
end

% North
j=Ny-1;
for i=2:Nx-2
        po=i+(j-1)*(Nx-1);
        A(po,po)= (-2/dx2)-(2/dy2); 
        A(po,po+1)=1/dx2; 
        A(po,po-1)=1/dx2;
        A(po,po-(Nx-1))=1/dy2;
    %    A(po,po+(Nx-1))=1/dy2;
end

% West
i=1;
for j=2:Ny-2
        po=i+(j-1)*(Nx-1);
        A(po,po)= (-2/dx2)-(2/dy2); 
        A(po,po+1)=1/dx2; 
     %   A(po,po-1)=1/dx2;
        A(po,po-(Nx-1))=1/dy2;
        A(po,po+(Nx-1))=1/dy2;
end

% East
i=Nx-1;
for j=2:Ny-2
        po=i+(j-1)*(Nx-1);
        A(po,po)= (-2/dx2)-(2/dy2); 
    %    A(po,po+1)=1/dx2; 
        A(po,po-1)=1/dx2;
        A(po,po-(Nx-1))=1/dy2;
        A(po,po+(Nx-1))=1/dy2;
end

% South-west
i=1;j=1;
        po=i+(j-1)*(Nx-1);
        A(po,po)= (-2/dx2)-(2/dy2); 
        A(po,po+1)=1/dx2; 
      %  A(po,po-1)=1/dx2;
       % A(po,po-(Nx-1))=1/dy2;
        A(po,po+(Nx-1))=1/dy2;
% South-east          
i=Nx-1;j=1;
        po=i+(j-1)*(Nx-1);
        A(po,po)= (-2/dx2)-(2/dy2); 
  %      A(po,po+1)=1/dx2; 
        A(po,po-1)=1/dx2;
   %     A(po,po-(Nx-1))=1/dy2;
        A(po,po+(Nx-1))=1/dy2;

% North-east         
i=Nx-1;j=Ny-1;
        po=i+(j-1)*(Nx-1);
        A(po,po)= (-2/dx2)-(2/dy2); 
      %  A(po,po+1)=1/dx2; 
        A(po,po-1)=1/dx2;
        A(po,po-(Nx-1))=1/dy2;
       % A(po,po+(Nx-1))=1/dy2;
% North-west          
i=1;j=Ny-1;
        po=i+(j-1)*(Nx-1);
        A(po,po)= (-2/dx2)-(2/dy2); 
        A(po,po+1)=1/dx2; 
 %       A(po,po-1)=1/dx2;
        A(po,po-(Nx-1))=1/dy2;
  %      A(po,po+(Nx-1))=1/dy2;

end