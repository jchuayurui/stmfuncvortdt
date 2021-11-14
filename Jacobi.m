function [u0]=Jacobi(A,b,u0,resobj)

% Jacobi method
D=diag(diag(A));
R=A-D;
k=0;

while 1
    u1=D\(b-R*u0); 
    res=norm(b-A*u1);
    if res<resobj
        break
    end
    u0=u1;
    k=k+1;
end

end
