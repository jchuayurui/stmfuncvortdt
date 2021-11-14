% A lid-driven cavity flow problem using Streamfunction-Vorticity formulation

clear;

% specify Reynolds number
Re=56/3;

% set up the grids
Nx=51;Ny=21;   % no. of intervals
Lx=3;Ly=1;

x=linspace(0,Lx,Nx+1);dx=x(2)-x(1);
y=linspace(0,Ly,Ny+1);dy=y(2)-y(1);
      % Here dx=dy, so we use dl to represent them

[x2,y2]=meshgrid(x,y);    % generate a 2D mesh

% specify the value of dt
%dt=0.015;   % dt_max = 0.00625
dt = 0.001
NMax=(2/dt);   % maximum of time steps. Total time = NMax*dt

% assemble A matrix to solve the Poisson equation Au=b
[A]=assembleA(Nx,Ny,dx,dy);

% initial condition for the vorticity
vort = zeros(Nx-1,Ny-1);
t=0;

%energyarray=[];
uarrayexp=[];
figure('units','normalized','outerposition',[0 0 1 1])

% time evoluion
for i=1:NMax

stmfunc = solve_Poisson(vort,A,Nx,Ny);            % solving the Poisson equation
                                                  % for streamfunction
vort = advance_vort(stmfunc,vort,Nx,Ny,dx,dy,dt,Re,t);
%vort = advance_vortimp(stmfunc,vort,Nx,Ny,dx,dy,dt,Re,t); % advance the vorticity transport 
                                                  % equation

disp(['Finish timestep' num2str(i)])

t=t+dt;

% for plotting
figure(1)
subplot(1,2,1)
[u,v]=get_uv(stmfunc,Nx,Ny,dx,dy,t);
quiver(x2,y2,u',v')
axis equal
xlim([0 3]);ylim([0 1.1]);
title(['At time step = ' num2str(i)])
set(gca,'FontSize',40)
xlabel('x');ylabel('y')
pause(0.2)
% drawnow

%energyarray = [energyarray sqrt(u(end-1,end-1)^2+v(end-1,end-1)^2)];
uarrayexp=[uarrayexp u(26,11)];
subplot(1,2,2)
plot(dt:dt:dt*i,uarrayexp,'-*b')
title(['Velocity u at (1.5,0.5)'])
set(gca,'FontSize',20)
xlabel('time');ylabel('u')

end
%vort = reshape(vort,Nx-1,Ny-1);

% post-processing
[u,v] = get_uv(stmfunc,Nx,Ny,dx,dy,t);

yexp = fft(uarrayexp);
fsexp = 1/dt;
fexp = (0:length(yexp)-1)*fsexp/length(yexp);
plot(fexp,abs(yexp));
xlabel('Frequency(Hz)')
ylabel('Magnitude')
title('Magnitude')
