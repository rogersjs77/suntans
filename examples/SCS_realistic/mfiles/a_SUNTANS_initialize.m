clear
% profile on

igrid = 0;
idepth = 1;
iprofile=0;
iinit_ic = 1;
iinit_bc = 1;

datadir='../rundata';
E = 401E3;
N =2120e3;
Zone = 50; % UTM zone
NCOM_file=['../../../../setup testing/from_Ko/NCOM_SCS_2300m/lzsnfs-201105.nc'];
BATHY_file = ['../../../../setup testing/SCS_bathymetry/bathy_1100m_4tyt.mat'];
% lon0 = 116.0110;
% lat0 = 20;
% [N,E,Zone,lcm]=ell2utm(deg2rad(lat0),deg2rad(lon0));
if igrid
    dx=2.3e3;
    L=475e3;
    W=390e3;%280e3;
    Nx = round(L/dx);%300;
    Ny = round(W/dx);%round(W/500);
    BC = [2 2 2 6]; %E,N,W,S Boundary conditions, 
    % 1 solid free-slip, 2 velocity, 3 free-surface, 5 periodic, 6 mixed 2/3
    STRETCHING=false;
    CHEBYCHEV=false;
    Lr=10;
    K=2;
    rmax=1.1;   
    x0=E;
    y0=N;
    theta=0;%-12;
    save SUNTANS_grid.mat -v7.3
%     profile on
    quadgrid_periodic(datadir,L,W,Nx,Ny,BC,STRETCHING,CHEBYCHEV,Lr,K,rmax,x0,y0,theta)
%     profile viewer
end

if idepth
    shape='NCOM';
    D0=25;
    x0 = [];
    y0 = [];
    theta=[];
    dtrend = [];
    bathy_dir = BATHY_file;
    dmin = 5000;
    extreme=0;
    SUNTANS_depth_ideal(datadir,D0,shape,bathy_dir,x0,y0,theta,dtrend,dmin,extreme)
end

if iprofile
    name = {'TC1_2011','M2_2013','M6_2013','S7_2005','B1_2005','EBC','NBC','WBC','SBC'};
    lat = [21.0707 21.01 20.83 21.614 21.365 ];
    lon = [117.2202 118.16 119.05 117.282 118.594 ];
    profile_points(datadir,name,lat,lon)
end

if iinit_ic
   initial_ideal = 0;
   init_filename = '/sun_IC.nc';
   SUNTANS_initial_conditions(initial_ideal,datadir,...
                 init_filename,NCOM_file);
end

if iinit_bc   
   bc_ideal = 0;
   bc_filename = '/sun_BC.nc';
   bc_dt = 3600; % ideal bc dt step [-1/24:1/24:4]*24*3600; % seconds since startime
   SUNTANS_boundary_conditions(bc_ideal,datadir,...
       bc_filename,bc_dt,NCOM_file);   
end
% profile viewer