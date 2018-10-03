function [] = SUNTANS_boundary_conditions(bc_ideal,datadir,...
                   bc_filename,bc_dt,NCOM_file)


%% %%%%%%%%
% SUNTANS Boundary Conditions
% Justin Rogers
% Stanford University
%% %%%%%%%

% clear
% bc_ideal = 1;
% datadir='../rundata';
% bcfilename = '/sun_BC.nc';
% bc_time = [0:100:1000]; % seconds since 1990-1-1
% UTM_zone = 50;
poolobj = parpool('local'); 
%% load grid data
clear mex

c = load([datadir,'/cells.dat']);
Nc = size(c,1);
if size(c,2)<9
    xv = c(:,1);
    yv = c(:,2);
else
    xv = c(:,2);
    yv = c(:,3);
end
cells = c(:,4:7);
neigh = c(:,8:11);
nfaces = c(:,1);
numsides =nfaces(1,1);
cellp=0:(Nc-1);

points = load([datadir,'/points.dat']);
Np = size(points,1);
xp = points(:,1);
yp = points(:,2);

ed = load([datadir,'/edges.dat']);
Ne = size(ed,1);
edges = ed(:,1:2);
mark = ed(:,3);
grad =ed(:,4:5);
edgep = 0:(Ne-1);
for j=1:Ne
    xe(j) = 0.5*(xp(edges(j,1)+1,1)+xp(edges(j,2)+1,1));
    ye(j) = 0.5*(yp(edges(j,1)+1,1)+yp(edges(j,2)+1,1));
    
    % get outward unit normal vectors for computational points only
    n1(j) = yp(edges(j,1)+1,1)-yp(edges(j,2)+1,1);
    n2(j) = xp(edges(j,1)+1,1)-xp(edges(j,2)+1,1);
    n = sqrt(n1(j).^2+n2(j).^2);
    n1(j) = n1(j)/n;
    n2(j) = n2(j)/n;
end
% we don't care about boundary points for sponge fluxes, set to 0
n1(mark>0)=0;
n2(mark>0)=0;

depth=load([datadir,'/depth.dat']);
dv = depth(:,3);

Nkmax = getvalue([datadir,'/suntans.dat'],'Nkmax');
Nk = Nkmax+zeros(Nc,1);
Nkw=Nkmax+1;

files = dir(['../data/vertspace.dat']);
if ~isempty(files) % load in suntans file
    dz = load(['../data/vertspace.dat']);
else % figure out myself
    rstretch = getvalue([datadir,'/suntans.dat'],'rstretch');
    dz=1;
    for i=1:Nkmax-1
        dz(i+1) = rstretch*dz(i);        
    end
    dz = dz'*max(dv)/sum(dz);
end
z_r = getz(dz); %depth

% get 2d normal vectors for easy computation
[~,N1]=meshgrid(z_r,n1);
[~,N2]=meshgrid(z_r,n2);

basetime = sprintf('%14.6f',getvalue([datadir,'/suntans.dat'],'basetime'));
mtime_base = datenum(basetime,'yyyymmdd.HHMMSS');
starttime = sprintf('%14.6f',getvalue([datadir,'/suntans.dat'],'starttime'));
mtime_start = datenum(starttime,'yyyymmdd');
Toffset = mtime_start-mtime_base; % days between basetime and starttime
dt = getvalue([datadir,'/suntans.dat'],'dt');
nsteps = getvalue([datadir,'/suntans.dat'],'nsteps');
mtime_end = mtime_start + dt*nsteps/86400; % ending time of simulation
ntaverage =  getvalue([datadir,'/suntans.dat'],'ntaverage');
ntaveragestore =  getvalue([datadir,'/suntans.dat'],'ntaveragestore');
sponge_distance =  getvalue([datadir,'/suntans.dat'],'sponge_distance');
sponge_distance =  getvalue([datadir,'/suntans.dat'],'sponge_distance');
Tfilt = dt*ntaverage*ntaveragestore/3600/24; % filter time window, days

beta =  getvalue([datadir,'/suntans.dat'],'beta');
gamma =  getvalue([datadir,'/suntans.dat'],'gamma');
rho0=1000;

%% convert to lat lon

GRID = load('./SUNTANS_grid');
if isfield(GRID,'Zone')
    UTM_zone = GRID.Zone;
else
    UTM_zone = 10;
end

[latv,lonv]=utm2ell(yv,xv,UTM_zone);
latv = rad2deg(latv);
lonv = rad2deg(lonv);

%% find indices of type2 and type3 bc

type2 = mark==2; % edges with type 2
type3_edge = grad(mark==3,1); % index of cells neighboring type 3 edge
type3 = logical(0*xv);
type3(type3_edge+1)=1;

Ntype2 = sum(type2);
Ntype3 = sum(type3);

%% get outward normal vector on edges
r_v = 0*xv';
rn1_v = 0*xv';
rn2_v = 0*xv';

for i=1:Nc    
  r2 = (xe-xv(i)).^2+(ye-yv(i)).^2;
  [r2,indx]=min(r2(type2));
  
  xb = xe(type2);
  xb = xb(indx);
  
  yb = ye(type2);
  yb = yb(indx);
  
  rn1_v(i) = (xb-xv(i))/sqrt(r2);
  rn2_v(i) = (yb-yv(i))/sqrt(r2);
  
  r_v(i) = sqrt((xb-xv(i)).^2+(yb-yv(i)).^2);

end

% interp to edges
F = scatteredInterpolant(xv,yv,rn1_v','linear','linear');
rn1_e = F(xe,ye);

F = scatteredInterpolant(xv,yv,rn2_v','linear','linear');
rn2_e = F(xe,ye);

F = scatteredInterpolant(xv,yv,r_v','linear','linear');
r_e = F(xe,ye);

% normalize so unit vector
r2 = sqrt(rn1_e.^2+rn2_e.^2);
if isempty(r2)
    disp('r2 is empty, not interpolating')
    return
end
rn1_e = rn1_e./r2;
rn2_e = rn2_e./r2;

D_hat = exp(-4.0*r_e/sponge_distance);

clear F r2

%%
disp('writing netcdf BC file')

% ensure no open/locked files
clear mex
ncid = netcdf.create([datadir bc_filename],'CLOBBER');
netcdf.close(ncid);
% delete(['..\rundata\sun_BC.nc'])

cmode = netcdf.getConstant('NETCDF4');
cmode = bitor(cmode,netcdf.getConstant('CLASSIC_MODEL'));
ncid = netcdf.create([datadir bc_filename],cmode);

netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Created', ['Created on ' datestr(now)]);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Author', '');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Description', 'SUNTANS Boundary Conditions File');

% define dimensions
Nk_dimid = netcdf.defDim(ncid,'Nk',Nkmax);
Nt_dimid = netcdf.defDim(ncid,'Nt',netcdf.getConstant('NC_UNLIMITED'));
Ne_dimid = netcdf.defDim(ncid,'Ne',Ne);
if Ntype2
    Ntype2_dimid = netcdf.defDim(ncid,'Ntype2',Ntype2);
end
if Ntype3
    Ntype3_dimid = netcdf.defDim(ncid,'Ntype3',Ntype3);
end

% define variables
varid = netcdf.defVar(ncid,'z','NC_DOUBLE',Nk_dimid);
netcdf.putAtt(ncid,varid,'long_name','Vertical grid mid-layer depth');
netcdf.putAtt(ncid,varid,'units','meters');

varid = netcdf.defVar(ncid,'time','NC_DOUBLE',Nt_dimid);
netcdf.putAtt(ncid,varid,'units','seconds since 1990-01-01 00:00:00');
netcdf.putAtt(ncid,varid,'long_name','Boundary time');
netcdf.defVarFill(ncid,varid,false,999999);

if Ntype2
    varid = netcdf.defVar(ncid,'xe','NC_DOUBLE',Ntype2_dimid);
    netcdf.putAtt(ncid,varid,'long_name','Easting of type-2 boundary points');
    netcdf.putAtt(ncid,varid,'units','meters');
    
    varid = netcdf.defVar(ncid,'ye','NC_DOUBLE',Ntype2_dimid);
    netcdf.putAtt(ncid,varid,'long_name','Northing of type-2 boundary points');
    netcdf.putAtt(ncid,varid,'units','meters');
    
    varid = netcdf.defVar(ncid,'xe_all','NC_DOUBLE',Ne_dimid);
    netcdf.putAtt(ncid,varid,'long_name','Easting of edge boundary points');
    netcdf.putAtt(ncid,varid,'units','meters');
    
    varid = netcdf.defVar(ncid,'ye_all','NC_DOUBLE',Ne_dimid);
    netcdf.putAtt(ncid,varid,'long_name','Northing of edge boundary points');
    netcdf.putAtt(ncid,varid,'units','meters');
    
    varid = netcdf.defVar(ncid,'edgep','NC_INT',Ntype2_dimid);
    netcdf.putAtt(ncid,varid,'long_name','Index of suntans grid edge corresponding to type-2 boundary');
    netcdf.putAtt(ncid,varid,'units','');
    
    varid = netcdf.defVar(ncid,'edgep_all','NC_INT',Ne_dimid);
    netcdf.putAtt(ncid,varid,'long_name','Index of suntans grid edge for all points');
    netcdf.putAtt(ncid,varid,'units','');
    
    varid = netcdf.defVar(ncid,'boundary_h','NC_DOUBLE',[Ntype2_dimid Nt_dimid]);
    netcdf.putAtt(ncid,varid,'long_name','Free-surface elevation at type-2 boundary point');
    netcdf.putAtt(ncid,varid,'units','meters');
    
    varid = netcdf.defVar(ncid,'boundary_u','NC_DOUBLE',[Ntype2_dimid Nk_dimid Nt_dimid]);
    netcdf.putAtt(ncid,varid,'long_name','Eastward velocity at type-2 boundary point');
    netcdf.putAtt(ncid,varid,'units','meters second-1');
    
    varid = netcdf.defVar(ncid,'boundary_v','NC_DOUBLE',[Ntype2_dimid Nk_dimid Nt_dimid]);
    netcdf.putAtt(ncid,varid,'long_name','Northward velocity at type-2 boundary point');
    netcdf.putAtt(ncid,varid,'units','meters second-1');
    
    varid = netcdf.defVar(ncid,'sponge_uf','NC_DOUBLE',[Ne_dimid Nk_dimid Nt_dimid]);
    netcdf.putAtt(ncid,varid,'long_name','Eastward sponge velocity at edges');
    netcdf.putAtt(ncid,varid,'units','meters second-1');
    
    varid = netcdf.defVar(ncid,'lowfreq_uf','NC_DOUBLE',[Ne_dimid Nk_dimid Nt_dimid]);
    netcdf.putAtt(ncid,varid,'long_name','Eastward low frequency velocity at edges');
    netcdf.putAtt(ncid,varid,'units','meters second-1');
    
    varid = netcdf.defVar(ncid,'sponge_u','NC_DOUBLE',[Ne_dimid Nk_dimid Nt_dimid]);
    netcdf.putAtt(ncid,varid,'long_name','Eastward sponge velocity at edges');
    netcdf.putAtt(ncid,varid,'units','meters second-1');
    
%     varid = netcdf.defVar(ncid,'sponge_v','NC_DOUBLE',[Ne_dimid Nk_dimid Nt_dimid]);
%     netcdf.putAtt(ncid,varid,'long_name','Northward sponge velocity at edges');
%     netcdf.putAtt(ncid,varid,'units','meters second-1');
%     
%     varid = netcdf.defVar(ncid,'lowfreq_u','NC_DOUBLE',[Ne_dimid Nk_dimid Nt_dimid]);
%     netcdf.putAtt(ncid,varid,'long_name','Eastward low frequency velocity at edges');
%     netcdf.putAtt(ncid,varid,'units','meters second-1');
%     
%     varid = netcdf.defVar(ncid,'lowfreq_v','NC_DOUBLE',[Ne_dimid Nk_dimid Nt_dimid]);
%     netcdf.putAtt(ncid,varid,'long_name','Northward low frequency velocity at edges');
%     netcdf.putAtt(ncid,varid,'units','meters second-1');
    
    varid = netcdf.defVar(ncid,'boundary_w','NC_DOUBLE',[Ntype2_dimid Nk_dimid Nt_dimid]);
    netcdf.putAtt(ncid,varid,'long_name','Vertical velocity at type-2 boundary point');
    netcdf.putAtt(ncid,varid,'units','meters second-1');
    
    varid = netcdf.defVar(ncid,'boundary_T','NC_DOUBLE',[Ntype2_dimid Nk_dimid Nt_dimid]);
    netcdf.putAtt(ncid,varid,'long_name','Water temperature at type-2 boundary point');
    netcdf.putAtt(ncid,varid,'units','degrees C');
    
    varid = netcdf.defVar(ncid,'boundary_S','NC_DOUBLE',[Ntype2_dimid Nk_dimid Nt_dimid]);
    netcdf.putAtt(ncid,varid,'long_name','Salinity at type-2 boundary point');
    netcdf.putAtt(ncid,varid,'units','psu');
end

if Ntype3
    varid = netcdf.defVar(ncid,'xv','NC_DOUBLE',Ntype3_dimid);
    netcdf.putAtt(ncid,varid,'long_name','Easting of type-3 boundary points');
    netcdf.putAtt(ncid,varid,'units','meters');
    
    varid = netcdf.defVar(ncid,'yv','NC_DOUBLE',Ntype3_dimid);
    netcdf.putAtt(ncid,varid,'long_name','Northing of type-3 boundary points');
    netcdf.putAtt(ncid,varid,'units','meters');
    
    varid = netcdf.defVar(ncid,'cellp','NC_INT',Ntype3_dimid);
    netcdf.putAtt(ncid,varid,'long_name','Index of suntans grid cell corresponding to type-3 boundary');
    netcdf.putAtt(ncid,varid,'units','');
    
    varid = netcdf.defVar(ncid,'uc','NC_DOUBLE',[Ntype3_dimid Nk_dimid Nt_dimid]);
    netcdf.putAtt(ncid,varid,'long_name','Eastward velocity at type-3 boundary point');
    netcdf.putAtt(ncid,varid,'units','meters second-1');
    
    varid = netcdf.defVar(ncid,'vc','NC_DOUBLE',[Ntype3_dimid Nk_dimid Nt_dimid]);
    netcdf.putAtt(ncid,varid,'long_name','Northward velocity at type-3 boundary point');
    netcdf.putAtt(ncid,varid,'units','meters second-1');
    
    varid = netcdf.defVar(ncid,'wc','NC_DOUBLE',[Ntype3_dimid Nk_dimid Nt_dimid]);
    netcdf.putAtt(ncid,varid,'long_name','Vertical velocity at type-3 boundary point');
    netcdf.putAtt(ncid,varid,'units','meters second-1');
    
    varid = netcdf.defVar(ncid,'T','NC_DOUBLE',[Ntype3_dimid Nk_dimid Nt_dimid]);
    netcdf.putAtt(ncid,varid,'long_name','Temperature at type-3 boundary point');
    netcdf.putAtt(ncid,varid,'units','degrees C');
    
    varid = netcdf.defVar(ncid,'S','NC_DOUBLE',[Ntype3_dimid Nk_dimid Nt_dimid]);
    netcdf.putAtt(ncid,varid,'long_name','Salinity at type-3 boundary point');
    netcdf.putAtt(ncid,varid,'units','psu');
    
    varid = netcdf.defVar(ncid,'h','NC_DOUBLE',[Ntype3_dimid Nt_dimid]);
    netcdf.putAtt(ncid,varid,'long_name','Water surface elevation at type-3 boundary point');
    netcdf.putAtt(ncid,varid,'units','meters');
    
   
end
% end define variables
netcdf.endDef(ncid);

% write grid variables
varid = netcdf.inqVarID(ncid,'z');
netcdf.putVar(ncid,varid,z_r);

if Ntype2
varid = netcdf.inqVarID(ncid,'xe');
netcdf.putVar(ncid,varid,xe(type2));

varid = netcdf.inqVarID(ncid,'ye');
netcdf.putVar(ncid,varid,ye(type2));

varid = netcdf.inqVarID(ncid,'xe_all');
netcdf.putVar(ncid,varid,xe);

varid = netcdf.inqVarID(ncid,'ye_all');
netcdf.putVar(ncid,varid,ye);

varid = netcdf.inqVarID(ncid,'edgep');
netcdf.putVar(ncid,varid,edgep(type2));

varid = netcdf.inqVarID(ncid,'edgep_all');
netcdf.putVar(ncid,varid,edgep);
end

if Ntype3
varid = netcdf.inqVarID(ncid,'xv');
netcdf.putVar(ncid,varid,xv(type3));

varid = netcdf.inqVarID(ncid,'yv');
netcdf.putVar(ncid,varid,yv(type3));

varid = netcdf.inqVarID(ncid,'cellp');
netcdf.putVar(ncid,varid,cellp(type3));
end

% close file
netcdf.close(ncid);


%% compute fields
disp('computing BC variables')

if bc_ideal %idealized initial conditions
    
    Nbc_end = ceil(dt*nsteps/bc_dt)+1; % add extra on end
    time = bc_dt*[-1:Nbc_end]+Toffset*86400;% time in s since basetime
    Nt = length(time);
    
    if Ntype2
        TM2 = getvalue([datadir '/suntans.dat'],'TM2')*3600; % convert to s
        ULm = getvalue([datadir '/suntans.dat'],'ULm');
        ULz = getvalue([datadir '/suntans.dat'],'ULz');
        UHtide = getvalue([datadir '/suntans.dat'],'UHtide');
        UHiw = getvalue([datadir '/suntans.dat'],'Uiw');
        VLm = getvalue([datadir '/suntans.dat'],'VLm');
        VLz = getvalue([datadir '/suntans.dat'],'VLz');
        Phitide = getvalue([datadir '/suntans.dat'],'Phitide');
        Phiiw = getvalue([datadir '/suntans.dat'],'Phiiw');
        Fhat_x = getvalue([datadir '/suntans.dat'],'Fhat_x');
        Fhat_y = getvalue([datadir '/suntans.dat'],'Fhat_y');
        alpha1 = getvalue([datadir '/suntans.dat'],'alpha1');
        alpha2 = getvalue([datadir '/suntans.dat'],'alpha2');
        delta = getvalue([datadir '/suntans.dat'],'delta');
        drho = getvalue([datadir '/suntans.dat'],'drho');
        h2 = getvalue([datadir '/suntans.dat'],'h2');
        D = getvalue([datadir '/suntans.dat'],'D');
        lambda = getvalue([datadir '/suntans.dat'],'lambda');
        
        F_hat = 0*xe; 
        F_hat(Fhat_x.*rn1_e + Fhat_y.*rn2_e<0)=1;
        [~,F_hat]=meshgrid(z_r,F_hat);
        
        kiw_x = 2*pi/lambda*Fhat_x;
        kiw_y = 2*pi/lambda*Fhat_y;
        kiw = sqrt(kiw_x^2+kiw_y^2);
        
        Ctide = sqrt(D*9.81);
        ktide = 2*pi/(TM2*Ctide);
        ktide_x = ktide*Fhat_x;
        ktide_y = ktide*Fhat_y;
        
        fz = alpha1*exp(-(-z_r+h2)/delta)+alpha2;
        [~,XE]=meshgrid(fz,xe);
        [FZ,YE]=meshgrid(fz,ye);
        
%         fx = exp(-(XE-max(xv)).^2/(2*0.2^2*(GRID.L/4)^2));
        % compute idealized velocities over entire domain          
        for k=1:Nt
            u_LS_bar(:,:,k) = ULm + ULz*FZ;
            U_LS_tilde(:,:,k) =  UHtide * (ktide_x/ktide) *...
                cos(ktide_x*XE + ktide_y*YE - 2*pi*time(k)/TM2 + Phitide);
            u_LS_tilde_prime(:,:,k) =  UHiw * FZ .* (kiw_x/kiw) .*...
                cos(kiw_x*XE + kiw_y*YE - 2*pi*time(k)/TM2 + Phiiw);
            
            v_LS_bar(:,:,k) = VLm + VLz*FZ;
            V_LS_tilde(:,:,k) =  UHtide * (ktide_y/ktide) *...
                cos(ktide_x*XE + ktide_y*YE - 2*pi*time(k)/TM2 + Phitide);
            v_LS_tilde_prime(:,:,k) =  UHiw * FZ .* (kiw_y/kiw) .*...
                cos(kiw_x*XE + kiw_y*YE - 2*pi*time(k)/TM2 + Phiiw);
        end       
        
        % velocity at boundary and sponges
        for k=1:Nt
            ubc(:,:,k) = u_LS_bar(:,:,k) + U_LS_tilde(:,:,k) + F_hat.*u_LS_tilde_prime(:,:,k);
            vbc(:,:,k) = v_LS_bar(:,:,k) + V_LS_tilde(:,:,k) + F_hat.*v_LS_tilde_prime(:,:,k);

            sponge_u(:,:,k) = F_hat.*u_LS_tilde_prime(:,:,k);
            sponge_v(:,:,k) = F_hat.*v_LS_tilde_prime(:,:,k);
            
            lowfreq_u(:,:,k) = u_LS_bar(:,:,k);
            lowfreq_v(:,:,k) = v_LS_bar(:,:,k);
            
            % face fluxes
            sponge_uf(:,:,k) = N1.*sponge_u(:,:,k) + N2.*sponge_v(:,:,k);
            lowfreq_uf(:,:,k) = N1.*lowfreq_u(:,:,k) + N2.*lowfreq_v(:,:,k);
        end
        
        
        % free surface
        boundary_h = zeros(Ntype2,Nt);

        % velocity
        boundary_u = ubc(type2,:,:);
        boundary_v = vbc(type2,:,:);
        boundary_w = 0*boundary_u;
        
        % salt
        boundary_S = 0*boundary_u;

        % temp
        a1 = 23.3637651923536;
        a2 = 3.18189680534772;
        a3 = -44.1196509264888;
        a4 = 293.118606564064;
        for i=1:Ntype2
            for j=1:Nt
                boundary_T(i,:,j) = a1*exp(-(-z_r+a3)/a4)+a2;
            end
        end
    end
    
    if Ntype3
        % eta
        h = zeros(Ntype3,Nt);
        
        % velocity
        uc = zeros(Ntype3,Nkmax,Nt);
        vc = uc;
        wc = uc;
        
        % salt
        S = uc;
        
        % temp
        a1 = 23.3637651923536;
        a2 = 3.18189680534772;
        a3 = -44.1196509264888;
        a4 = 293.118606564064;
        for i=1:Ntype3
            for j=1:Nt
                T(i,:,j) = a1*exp(-(-z_r+a3)/a4)+a2;
            end
        end   
        
    end

else % realistic initial conditions from NCOM
    
    ncid = netcdf.open(NCOM_file,'NC_NOWRITE');
    fill_value = -99;
    
    varid = netcdf.inqDimID(ncid,'time');
    [~, tcount] = netcdf.inqDim(ncid,varid);
    tstart=0; % starts at 0
    
    % time
    varname = netcdf.inqVar(ncid,0);
    varid = netcdf.inqVarID(ncid,varname);
    NC_time_all = netcdf.getVar(ncid,varid,tstart,tcount,'double');
    NCOM_basetime = datenum(netcdf.getAtt(ncid,varid,'start'));
    time_step = netcdf.getAtt(ncid,varid,'units');
    if strcmp(time_step,'hour')
        time_step = 24;
    end
    NC_mtime_all = NC_time_all/time_step+NCOM_basetime;

    % lon
    varname = netcdf.inqVar(ncid,1);
    varid = netcdf.inqVarID(ncid,varname);
    NC_lon = netcdf.getVar(ncid,varid,'double');

    % lat
    varname = netcdf.inqVar(ncid,2);
    varid = netcdf.inqVarID(ncid,varname);
    NC_lat = netcdf.getVar(ncid,varid,'double');

    %depth
    varname = netcdf.inqVar(ncid,3);
    varid = netcdf.inqVarID(ncid,varname);
    NC_depth = netcdf.getVar(ncid,varid,'double');
    NC_z = -NC_depth;

    NC_xcount = length(NC_lon);
    NC_ycount = length(NC_lat);
    NC_zcount = length(NC_depth);
    
    [NC_LAT,NC_LON]=meshgrid(NC_lat,NC_lon);    
    [NC_Y,NC_X,Zone,lcm]=ell2utm(deg2rad(NC_LAT),deg2rad(NC_LON));
    lcm = nanmedian(reshape(lcm(Zone==GRID.Zone),[],1));
    [NC_Y,NC_X,~,~]=ell2utm(deg2rad(NC_LAT),deg2rad(NC_LON),lcm);
    
    NC_dy = max(max(abs(diff(NC_Y,[],2))));
    NC_dx = max(max(abs(diff(NC_X,[],1))));
    
    % get buffer points, ccw from ll, dx from model 
    buffer =round(sponge_distance/NC_dx);
    x1 = min(xe)-buffer*NC_dx;
    x2 = max(xe)+buffer*NC_dx;    
    y1 = min(ye)-buffer*NC_dy;
    y2 = max(ye)+buffer*NC_dy;    
    xbuffer = [x1 x2 x2 x1 x1];
    ybuffer = [y1 y1 y2 y2 y1];
    
    % find spatial indices (NC files start at 0)
    indx = inpolygon(NC_X,NC_Y,xbuffer,ybuffer);
    % rows are x, columns are y
    xsum = find(sum(indx,2)>0);
    xstart=xsum(1)-1;
    xend =xsum(end)-1;
    xcount = xend-xstart+1;
    % rows are x, columns are y
    ysum = find(sum(indx,1)>0);
    ystart=ysum(1)-1;
    yend =ysum(end)-1;
    ycount = yend-ystart+1;
    
    % find z index
    zcount = length(NC_depth);
    
    % find time index closest to starttime
    [~,tstart] = min(abs(mtime_start-NC_mtime_all));
    if NC_mtime_all(tstart)>=mtime_start % ensure start indx is before starttime
       tstart=tstart-1; 
    end
    tstart = tstart-1; %netcdf indexing
    % find time index closest to endttime
    [~,tend] = min(abs(mtime_end-NC_mtime_all));
    if NC_mtime_all(tend)<=mtime_end % ensure end indx is after starttime
       tend=tend+1; 
    end
    tend = tend+1; %netcdf indexing, add 1, so 2 nc points after last step
    tcount = tend-tstart+1;
    % find low pass filter tcount, closest to bin window
    tcount_filt = round(Tfilt/median(diff(NC_mtime_all)));
    
    % reload time, lat, lon data only in area of interest
    
    % time
    varid = netcdf.inqVarID(ncid,'time');
    NC_time = netcdf.getVar(ncid,varid,tstart,tcount,'double');
    NC_mtime =  NC_time/24+NCOM_basetime;
    time = (NC_mtime-mtime_base)*86400; % time in s since basetime
    
    % lon
    varname = netcdf.inqVar(ncid,1);
    varid = netcdf.inqVarID(ncid,varname);
    NC_lon = netcdf.getVar(ncid,varid,[xstart],[xcount],'double');

    % lat
    varname = netcdf.inqVar(ncid,2);
    varid = netcdf.inqVarID(ncid,varname);
    NC_lat = netcdf.getVar(ncid,varid,ystart,ycount,'double');
    
    [NC_LAT,NC_LON]=meshgrid(NC_lat,NC_lon);    
    [NC_Y,NC_X,~,~]=ell2utm(deg2rad(NC_LAT),deg2rad(NC_LON),lcm);
    
    % get outward normal vectors, interp from model grid
    F = scatteredInterpolant(xv,yv,rn1_v','linear','linear');
    NC_rn1 = F(NC_X,NC_Y);
    F = scatteredInterpolant(xv,yv,rn2_v','linear','linear');
    NC_rn2 = F(NC_X,NC_Y);
    % normalize so unit vector
    r2 = sqrt(NC_rn1.^2+NC_rn2.^2);
    NC_rn1 = NC_rn1./r2;
    NC_rn2 = NC_rn2./r2;
    clear F r2
    
    % compute dz, z on faces
    NC_zface = [0; 0.5*(NC_z(1:end-1)+NC_z(2:end))];
    NC_zface = [NC_zface; NC_z(end)+(NC_z(end)-NC_zface(end))];

    NC_dz = diff(NC_zface);
    % make 3D dZ for easy computation
    for i=1:length(NC_lon)
        for j=1:length(NC_lat)
            for k=1:length(NC_z)
               NC_dZ(i,j,k)=NC_dz(k); 
            end
        end
    end
    grav=9.81;
    
    tindxi = tstart; %current time index
    tindx_lp = tindxi-tcount_filt;
    tcount_load = tcount_filt+1+tcount; % time steps to load        
    if tindx_lp<0 % check that time data exists
        tindx_lp=0;
        tcount_load = tindxi+1+tcount;
        disp('warning, not enough data available for full time filter')
    end    
    
    % load in netcdf data within space and time limits for analysis
    varid = netcdf.inqVarID(ncid,'time');
    NC_time_lp = netcdf.getVar(ncid,varid,tindx_lp,tcount_load,'double');
    NC_mtime_lp =  NC_time_lp/24+NCOM_basetime;
    
    disp('reading netcdf BC data')
    varid = netcdf.inqVarID(ncid,'ssh');
    NC_ssh = netcdf.getVar(ncid,varid,[xstart ystart tindx_lp],[xcount ycount tcount_load],'double');
    NC_ssh(NC_ssh==fill_value)=nan;

    varid = netcdf.inqVarID(ncid,'temperature');
    NC_temp = netcdf.getVar(ncid,varid,[xstart ystart 0 tindx_lp],[xcount ycount NC_zcount tcount_load],'double');
    NC_temp(NC_temp==fill_value)=nan;

    varid = netcdf.inqVarID(ncid,'salinity');
    NC_salt = netcdf.getVar(ncid,varid,[xstart ystart 0 tindx_lp],[xcount ycount NC_zcount tcount_load],'double');
    NC_salt(NC_salt==fill_value)=nan;

    varid = netcdf.inqVarID(ncid,'u-velocity');
    NC_u = netcdf.getVar(ncid,varid,[xstart ystart 0 tindx_lp],[xcount ycount NC_zcount tcount_load],'double');
    NC_u(NC_u==fill_value)=nan;

    varid = netcdf.inqVarID(ncid,'v-velocity');
    NC_v = netcdf.getVar(ncid,varid,[xstart ystart 0 tindx_lp],[xcount ycount NC_zcount tcount_load],'double');
    NC_v(NC_v==fill_value)=nan;
    
    % compute density
    NC_rho = rho0*(1+beta*NC_salt-gamma*NC_temp);
    
    % temporary values for parallel interp
    xe_type2 = xe(type2);
    ye_type2 = ye(type2);
    xv_type3 = xv(type3);
    yv_type3 = yv(type3);
    
    

    for i=1:tcount
        disp(['interpolating BC from NCOM netcdf ' datestr(NC_mtime(i))])
        
        clear p_prime u v rho_prime U V u_prime v_prime U_tilde V_tilde Fx_prime Fy_prime
        % find index of current step on lp time series
        [~,tindxi] = min(abs(NC_mtime(i)-NC_mtime_lp));
        % index for averaging
        tindx_avg = (tindxi-tcount_filt):1:tindxi;
        tindx_avg = tindx_avg(tindx_avg>0);
        
        % ssh
        ssh = NC_ssh(:,:,tindxi);
        %% decompose free surface
        % spatial linear fit to surface is tidal component
        ydata = reshape(ssh,[],1);
        xdata = reshape(NC_X,[],1);
        xdata2 = reshape(NC_Y,[],1);
        indxe = ~isnan(ydata); % in case any dry cells
        [px]=polyfit(xdata(indxe),ydata(indxe),1);
        [py]=polyfit(xdata2(indxe),ydata(indxe),1);
        etatide = px(1)*NC_X + py(1)*NC_Y + px(2);
        etawave = ssh - etatide;
        etawave = etawave - nanmean(nanmean(etawave));
        % remove the wave component for stability
        ssh_smooth = ssh-etawave; 
        
        % temp
        temp = NC_temp(:,:,:,tindxi);
        
        % salt
        salt = NC_salt(:,:,:,tindxi);
        
        tic
        % background density
        rho_b = squeeze(nanmean(NC_rho(:,:,:,tindx_avg),4));
        for k=1:length(tindx_avg)
           rho_prime(:,:,:,k)=NC_rho(:,:,:,tindx_avg(k))-rho_b;            
        end
        
        % pressure
        p_prime = zeros(length(NC_lon),length(NC_lat),length(NC_depth),length(tindx_avg));
        for k=1:length(tindx_avg)
            p_prime(:,:,1,k) = 0.5*squeeze(rho_prime(:,:,1,k)).*grav.*(NC_dz(1)+NC_ssh(:,:,tindx_avg(k))); % including the free surface in the pressure
        end
        %    here's a p_prime calculation using cut cell dz values at the bottom
        for k = 2:length(NC_depth)
            p_prime(:,:,k,:) = p_prime(:,:,k-1,:)+0.5*rho_prime(:,:,k-1,:).*grav.*NC_dZ(:,:,k-1)+...
                               0.5*rho_prime(:,:,k,:).*grav.*NC_dZ(:,:,k);
        end
        
        % velocity decomposition
        u = NC_u(:,:,:,tindx_avg);
        v = NC_v(:,:,:,tindx_avg);
        
        % barotropic
        for k=1:length(tindx_avg)
            U(:,:,k) = nansum(squeeze(u(:,:,:,k)).*NC_dZ,3)./nansum(NC_dZ,3);
            V(:,:,k) = nansum(squeeze(v(:,:,:,k)).*NC_dZ,3)./nansum(NC_dZ,3);
        end        

        % baroclinic
        for k=1:size(u,3)
            u_prime(:,:,k,:) = squeeze(u(:,:,k,:))-U;
            v_prime(:,:,k,:) = squeeze(v(:,:,k,:))-V;
        end
        
        % mean barotropic
        Umean = nanmean(U,3);
        Vmean = nanmean(V,3);

        % barotropic high (tide)
        for k=1:length(tindx_avg)
           U_tilde(:,:,k)=U(:,:,k)-Umean;
           V_tilde(:,:,k)=V(:,:,k)-Vmean;
        end

        % mean baroclinic
        u_mean_prime = nanmean(u_prime,4);
        v_mean_prime = nanmean(v_prime,4);

        % baroclinic high (wave)
        for k=1:length(tindx_avg)
            u_tilde_prime(:,:,:,k) = u_prime(:,:,:,k)-u_mean_prime;
            v_tilde_prime(:,:,:,k) = v_prime(:,:,:,k)-v_mean_prime;
        end
        
        % estimate wave energy flux
        for k=1:length(tindx_avg)
            Fx_prime(:,:,k) = nansum(u_prime(:,:,:,k).*p_prime(:,:,:,k).*NC_dZ,3);
            Fy_prime(:,:,k) = nansum(v_prime(:,:,:,k).*p_prime(:,:,:,k).*NC_dZ,3);
        end
        Fx_prime_mean = nanmean(Fx_prime,3);
        Fy_prime_mean = nanmean(Fy_prime,3);
        
        % downsample to size of sponge D to eliminate small islands and smooth transition
        DS_scale = 2*2*floor(sponge_distance/NC_dx/2)+1; % ensure odd        
%         [DS_X,DS_Y,DS_Fx]=downsample_grid(NC_X,NC_Y,Fx_prime_mean,DS_scale);
%         [DS_X,DS_Y,DS_Fy]=downsample_grid(NC_X,NC_Y,Fy_prime_mean,DS_scale);
%         
        % interpolate flux back to NC grid
%         F = scatteredInterpolant(reshape(DS_X,[],1),reshape(DS_Y,[],1),reshape(DS_Fx,[],1));
%         Fx_lp = F(NC_X,NC_Y);
%         F = scatteredInterpolant(reshape(DS_X,[],1),reshape(DS_Y,[],1),reshape(DS_Fy,[],1));
%         Fy_lp = F(NC_X,NC_Y);
%         clear F
        
        hfilt = ones(DS_scale,DS_scale)/DS_scale^2;
        Fx_lp = imfilter(Fx_prime_mean,hfilt);
        Fy_lp = imfilter(Fy_prime_mean,hfilt);
        F_lp = sqrt(Fx_lp.^2+Fy_lp.^2);
        
        % set small values to zero
%         Fthresh = quantile(reshape(F_lp,[],1),0.5);
%         indx = F_lp<Fthresh;
%         Fx_lp(indx)=0;
%         Fy_lp(indx)=0;
        
        % compute Fhat
        % Fhat = [F dot n < 0] = 1, else 0
        Fhat = 0*Fx_prime_mean;
        Fhat(Fx_lp.*NC_rn1 + Fy_lp.*NC_rn2<0)=1;
        Fhat_lp = Fhat;
        
        toc
%         Fhat = 0*Fx_prime_mean;
%         Fhat(Fx.*NC_rn1 + Fy.*NC_rn2<0)=1;
%         Fhat_lp2 = imfilter(Fhat,hfilt);
%         Fhat_lp2(Fhat_lp2>0.5)=1;
%         Fhat_lp2(Fhat_lp2<=0.5)=0;
        
%         Fhat(Fx_prime_mean.*NC_rn1 + Fy_prime_mean.*NC_rn2<0)=1;  
%         [DS_X,DS_Y,DS_Fhat]=downsample_grid(NC_X,NC_Y,Fhat,DS_scale);
        
        % interpolate flux back to NC grid
%         F = scatteredInterpolant(reshape(DS_X,[],1),reshape(DS_Y,[],1),reshape(DS_Fhat,[],1));
%         Fhat_lp = F(NC_X,NC_Y);
        % now compute variables for writing to bc files
        % boundary velocities
        for k=1:length(NC_z)
            ubc(:,:,k) = U(:,:,end) + u_mean_prime(:,:,k) + Fhat_lp.*u_tilde_prime(:,:,k,end);
            vbc(:,:,k) = V(:,:,end) + v_mean_prime(:,:,k) + Fhat_lp.*v_tilde_prime(:,:,k,end);
            
            ulowfreq(:,:,k) = Umean + u_mean_prime(:,:,k);
            vlowfreq(:,:,k) = Vmean + v_mean_prime(:,:,k);
            
        end
        
        % sponge velocities
        for k=1:length(NC_z)
           usponge(:,:,k) = Fhat_lp.*u_tilde_prime(:,:,k,end);
           vsponge(:,:,k) = Fhat_lp.*v_tilde_prime(:,:,k,end);
        end
        
        %% now use parallel pool for faster interpolation

%         disp('interpolating NCOM variables to SUNTANS grid')
        for j=1
            % now interp variables to model grid 
            if Ntype2
                fprintf('Type 2: eta...')            
                boundary_h(:,i) = interpNCOM_parallel(NC_X,NC_Y,NC_z,ssh_smooth,xe_type2,ye_type2,z_r);
                
                fprintf('u...')
                boundary_u(:,:,i) = interpNCOM_parallel(NC_X,NC_Y,NC_z,ubc,xe_type2,ye_type2,z_r);
                
                fprintf('v...') 
                boundary_v(:,:,i) = interpNCOM_parallel(NC_X,NC_Y,NC_z,vbc,xe_type2,ye_type2,z_r);
                
                boundary_w(:,:,i) = 0*boundary_u(:,:,i);
                
                fprintf('sponge_u...')
                sponge_u(:,:,i) = interpNCOM_parallel(NC_X,NC_Y,NC_z,usponge,xe,ye,z_r);
                
                fprintf('sponge_v...') 
                sponge_v(:,:,i) = interpNCOM_parallel(NC_X,NC_Y,NC_z,vsponge,xe,ye,z_r);
                
                fprintf('lowfreq_u...')
                lowfreq_u(:,:,i) = interpNCOM_parallel(NC_X,NC_Y,NC_z,ulowfreq,xe,ye,z_r);
                
                fprintf('lowfreq_v...') 
                lowfreq_v(:,:,i) = interpNCOM_parallel(NC_X,NC_Y,NC_z,vlowfreq,xe,ye,z_r);
                
                % compute face velocities
                sponge_uf(:,:,i) = N1.*sponge_u(:,:,i) + N2.*sponge_v(:,:,i);
                lowfreq_uf(:,:,i) = N1.*lowfreq_u(:,:,i) + N2.*lowfreq_v(:,:,i);
                
                fprintf('T...')
                boundary_T(:,:,i) = interpNCOM_parallel(NC_X,NC_Y,NC_z,temp,xe_type2,ye_type2,z_r);
                
                fprintf('S...\n') 
                boundary_S(:,:,i) = interpNCOM_parallel(NC_X,NC_Y,NC_z,salt,xe_type2,ye_type2,z_r);
            end
            if Ntype3
                fprintf('Type 3: eta...')
                h(:,i) = interpNCOM_parallel(NC_X,NC_Y,NC_z,ssh_smooth,xv_type3,yv_type3,z_r);
                fprintf('uc...')
                uc(:,:,i) = interpNCOM_parallel(NC_X,NC_Y,NC_z,ubc,xv_type3,yv_type3,z_r);
                fprintf('vc...')
                vc(:,:,i) = interpNCOM_parallel(NC_X,NC_Y,NC_z,vbc,xv_type3,yv_type3,z_r);
                wc(:,:,i) = 0*uc(:,:,i);
%                 sponge_uc(:,:,i) = interpNCOM_parallel(NC_X,NC_Y,NC_z,usponge,xv,yv,z_r);
%                 sponge_vc(:,:,i) = interpNCOM_parallel(NC_X,NC_Y,NC_z,vsponge,xv,yv,z_r);
%                 lowfreq_v(:,:,i) = interpNCOM_parallel(NC_X,NC_Y,NC_z,vlowfreq,xe,ye,z_r);
                fprintf('T...')
                T(:,:,i) = interpNCOM_parallel(NC_X,NC_Y,NC_z,temp,xv_type3,yv_type3,z_r); 
                fprintf('S...')
                S(:,:,i) = interpNCOM_parallel(NC_X,NC_Y,NC_z,salt,xv_type3,yv_type3,z_r);
            end
                                    
        end       
   
        fprintf('done!\n')
 
        %% now assign values to BC arrays
        
        BC.Fhat(:,:,i) = Fhat_lp;
        BC.Fx(:,:,i) = Fx_prime_mean;
        BC.Fy(:,:,i) = Fy_prime_mean;
        BC.eta(:,:,i) = NC_ssh(:,:,tindxi);
        
    end   
    
    BC.X = NC_X;
    BC.Y = NC_Y;
    BC.xv = xv;
    BC.yv = yv;
    BC.mtime = NC_mtime;
    save('./SUNTANS_BC.mat','BC')
    
end

%% write variables
disp('writing variables to netcdf BC file')

ncid = netcdf.open([datadir bc_filename],'WRITE');

varid = netcdf.inqVarID(ncid,'time');
netcdf.putVar(ncid, varid,0,length(time), time);

if Ntype2
    varid = netcdf.inqVarID(ncid,'boundary_h');
    netcdf.putVar(ncid,varid,boundary_h);
    
    varid = netcdf.inqVarID(ncid,'boundary_u');
    netcdf.putVar(ncid,varid,boundary_u);
    
    varid = netcdf.inqVarID(ncid,'boundary_v');
    netcdf.putVar(ncid,varid,boundary_v);
    
    varid = netcdf.inqVarID(ncid,'sponge_uf');
    netcdf.putVar(ncid,varid,sponge_uf);
    
    varid = netcdf.inqVarID(ncid,'lowfreq_uf');
    netcdf.putVar(ncid,varid,lowfreq_uf);
    
%     varid = netcdf.inqVarID(ncid,'sponge_u');
%     netcdf.putVar(ncid,varid,sponge_u);
%     
%     varid = netcdf.inqVarID(ncid,'sponge_v');
%     netcdf.putVar(ncid,varid,sponge_v);
%     
%     varid = netcdf.inqVarID(ncid,'lowfreq_u');
%     netcdf.putVar(ncid,varid,lowfreq_u);
%     
%     varid = netcdf.inqVarID(ncid,'lowfreq_v');
%     netcdf.putVar(ncid,varid,lowfreq_v);
    
    varid = netcdf.inqVarID(ncid,'boundary_w');
    netcdf.putVar(ncid,varid,boundary_w);
    
    varid = netcdf.inqVarID(ncid,'boundary_T');
    netcdf.putVar(ncid,varid,boundary_T);
    
    varid = netcdf.inqVarID(ncid,'boundary_S');
    netcdf.putVar(ncid,varid,boundary_S);    
end

if Ntype3 
    varid = netcdf.inqVarID(ncid,'h');
    netcdf.putVar(ncid,varid,h);
    
    varid = netcdf.inqVarID(ncid,'uc');
    netcdf.putVar(ncid,varid,uc);
    
    varid = netcdf.inqVarID(ncid,'vc');
    netcdf.putVar(ncid,varid,vc);
    
    varid = netcdf.inqVarID(ncid,'wc');
    netcdf.putVar(ncid,varid,wc);
    
    varid = netcdf.inqVarID(ncid,'T');
    netcdf.putVar(ncid,varid,T);
    
    varid = netcdf.inqVarID(ncid,'S');
    netcdf.putVar(ncid,varid,S);    
end

netcdf.close(ncid);

%% plot bc
if bc_ideal
   subplot(3,1,1)
   uplot = squeeze(sqrt(ubc(:,1,end).^2+vbc(:,1,end).^2));
   scatter(xe,ye,5,uplot)
   cb=colorbar;
   ylabel(cb,'|ubc| top')
   
   subplot(3,1,2)
   uplot = squeeze(sqrt(sponge_u(:,1,end).^2+sponge_v(:,1,end).^2));
   scatter(xe,ye,5,uplot)
   cb=colorbar;
   ylabel(cb,'|usponge| top')
   ylabel('y (m)') 
   
   
   subplot(3,1,3)
   uplot = squeeze(sqrt(lowfreq_u(:,1,end).^2+lowfreq_v(:,1,end).^2));
   scatter(xe,ye,5,uplot)
   cb=colorbar;
   ylabel(cb,'|u LF| top')
   ylabel('y (m)')
   xlabel('x (m)')  
   
   
   print -djpeg -r300 figure_BC
    
    close
end


if bc_ideal==0
    
    qx = ceil(size(BC.Fx,2)/15);
    
    figure(345)
    subplot(3,1,1)
    title('eta with wave flux F'', last time step')
    pcolorjw(BC.X/1000,BC.Y/1000,sqrt(BC.Fx(:,:,end).^2+BC.Fy(:,:,end).^2))
    hold on
    
    f1 = BC.Fx(:,:,end);
    f2 = BC.Fy(:,:,end);
    quiver(BC.X(1:qx:end,1:qx:end)/1000,BC.Y(1:qx:end,1:qx:end)/1000,f1(1:qx:end,1:qx:end),f2(1:qx:end,1:qx:end),'color','k')
    cb=colorbar;
    ylabel(cb,'F''')
    colormap(dark_french)
    ylabel('y UTM (m)')
    axis equal
    
    subplot(3,1,2)
    pcolorjw(BC.X/1000,BC.Y/1000,BC.Fhat(:,:,end))
    cb= colorbar;
    ylabel(cb,'Fhat')
    ylabel('y UTM (m)')
    xlabel('x UTM (km)')
    axis equal
    
    subplot(3,1,3)
    scatter(xe/1000,ye/1000,5,D_hat)
    axis tight
    cb=colorbar;
    caxis([0 1])
    ylabel(cb,'Dhat')
    ylabel('y UTM (m)')
    xlabel('x UTM (km)')
    axis equal
    
    print -djpeg -r300 figure_BC
    
    close
%     cc = redbluecmap(64);
%     for i=1:3
%     ccc(:,i) = interp1(linspace(1,64,11),cc(:,i),1:64);
%     end
%     colormap(ccc)
    
    
end

delete(poolobj);