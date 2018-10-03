clear

igrid = 1;
idepth = 1;

datadir='../rundata';
if igrid
    dx=1000;
    L=4*130e3;
    W=3*dx;
    Nx = round(L/dx);%300;
    Ny = 3;%round(W/500);
    BC = [2 5 2 5]; %E,N,W,S Boundary conditions, 
    % 1 solid free-slip, 2 velocity, 3 free-surface, 5 periodic
    STRETCHING=false;
    CHEBYCHEV=false;
    Lr=10;
    K=2;
    rmax=1.1;
    save SUNTANS_grid.mat
 
    quadgrid_periodic(datadir,L,W,Nx,Ny,BC,STRETCHING,CHEBYCHEV,Lr,K,rmax)
    
end

if idepth
    shape='const';
    D0=3000;
    x0 = 132;
    y0 = 84;
    theta=56;
    dtrend = 'yes';
    bathy_dir = [];
    SUNTANS_depth_ideal(datadir,D0,shape,bathy_dir,x0,y0,theta,dtrend)
end