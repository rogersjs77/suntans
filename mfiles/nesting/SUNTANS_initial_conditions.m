function [] = SUNTANS_initial_conditions(initial_ideal,datadir,...
                 init_filename,NCOM_file,etainterp)
%%
% SUNTANS Initial Conditions
% Justin Rogers
% Stanford Univerity
%% %%%



% initial_ideal = 1;
% datadir='../rundata';
% initfilename = '/sun_IC.nc';
% init_time = 0; % seconds since 1990-1-1
% UTM_zone = 50;
delete(gcp('nocreate'))
poolobj = parpool('local'); 
%% load grid data

if nargin<5
    etainterp=0; % let eta_ic=0
end

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

points = load([datadir,'/points.dat']);
Np = size(points,1);
xp = points(:,1);
yp = points(:,2);

ed = load([datadir,'/edges.dat']);
Ne = size(ed,1);
edges = ed(:,1:2);
mark = ed(:,3);
grad =ed(:,4:5);
edgep = 0:(Ne-1); % this is the edge index in C format
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
n2=-n2; % set correct direction
% we don't care about boundary points for sponge fluxes, set to 0
n1(mark>0)=0;
n2(mark>0)=0;


depth=load([datadir,'/depth.dat']);
dv = depth(:,3);

Nkmax = getvalue([datadir,'/suntans.dat'],'Nkmax');
Nk = Nkmax+zeros(Nc,1);
Nkw=Nkmax+1;
basetime = sprintf('%14.6f',getvalue([datadir,'/suntans.dat'],'basetime'));
mtime_base = datenum(basetime,'yyyymmdd.HHMMSS');
starttime = sprintf('%14.6f',getvalue([datadir,'/suntans.dat'],'starttime'));
mtime_start = datenum(starttime,'yyyymmdd.HHMMSS');
Toffset = mtime_start-mtime_base; % days between basetime and starttime

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


save('SUNTANS_grid.mat','z_r','dz','-append')

% time = init_time;

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

%% create netcdf file
disp('creating netcdf IC file');

clear mex
delete([datadir '/*sun_IC.nc']);

cmode = netcdf.getConstant('NETCDF4');
cmode = bitor(cmode,netcdf.getConstant('CLASSIC_MODEL'));
ncid = netcdf.create([datadir init_filename],cmode);

netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Created', ['Created on ' datestr(now)]);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Author', '');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Description', 'SUNTANS Initial Conditions File');

% define dimensions
Nc_dimid = netcdf.defDim(ncid,'Nc',Nc);
Np_dimid = netcdf.defDim(ncid,'Np',Np);
Ne_dimid = netcdf.defDim(ncid,'Ne',Ne);
Nk_dimid = netcdf.defDim(ncid,'Nk',Nkmax);
Nkw_dimid = netcdf.defDim(ncid,'Nkw',Nkw);
numsides_dimid = netcdf.defDim(ncid,'numsides',numsides);
two_dimid = netcdf.defDim(ncid,'two',2);
time_dimid = netcdf.defDim(ncid,'time',netcdf.getConstant('NC_UNLIMITED'));

% define variables
varid = netcdf.defVar(ncid,'Nk','NC_INT',Nc_dimid);
netcdf.putAtt(ncid,varid,'long_name','Number of layers at face');

varid = netcdf.defVar(ncid,'neigh','NC_INT',[numsides_dimid Nc_dimid ]);
netcdf.putAtt(ncid,varid,'long_name','Maps every face to its neighbouring faces');
netcdf.putAtt(ncid,varid,'cf_role','face_face_connectivity');

varid = netcdf.defVar(ncid,'xp','NC_DOUBLE',Np_dimid);
netcdf.putAtt(ncid,varid,'long_name','Easting of 2D mesh node');
netcdf.putAtt(ncid,varid,'standard_name','Easting');

varid = netcdf.defVar(ncid,'xv','NC_DOUBLE',Nc_dimid);
netcdf.putAtt(ncid,varid,'long_name','Easting of 2D mesh face');
netcdf.putAtt(ncid,varid,'standard_name','Easting');

varid = netcdf.defVar(ncid,'z_r','NC_DOUBLE',Nk_dimid);
netcdf.putAtt(ncid,varid,'long_name','depth at layer mid points');
netcdf.putAtt(ncid,varid,'standard_name','ocean_z_coordinate');
netcdf.putAtt(ncid,varid,'units','m');
netcdf.putAtt(ncid,varid,'positive','up');

varid = netcdf.defVar(ncid,'latv','NC_DOUBLE',Nc_dimid);
netcdf.putAtt(ncid,varid,'long_name','Latitude of 2D mesh face');
netcdf.putAtt(ncid,varid,'standard_name','Latitude');

varid = netcdf.defVar(ncid,'mark','NC_INT',Ne_dimid);
netcdf.putAtt(ncid,varid,'long_name','Edge marker type');
netcdf.putAtt(ncid,varid,'units','0 - comp, 1 - boundary, 2, 3');
netcdf.putAtt(ncid,varid,'coordinates','xe, ye');

varid = netcdf.defVar(ncid,'lonv','NC_DOUBLE',Nc_dimid);
netcdf.putAtt(ncid,varid,'long_name','Longitude of 2D mesh face');
netcdf.putAtt(ncid,varid,'standard_name','Longitude');

varid = netcdf.defVar(ncid,'edges','NC_INT',[two_dimid Ne_dimid]);
netcdf.putAtt(ncid,varid,'long_name','Maps every edge to the two nodes it connects');
netcdf.putAtt(ncid,varid,'cf_role','edge_node_connectivity');

varid = netcdf.defVar(ncid,'dz','NC_DOUBLE',Nk_dimid);
netcdf.putAtt(ncid,varid,'long_name','z layer spacing');
netcdf.putAtt(ncid,varid,'units','m');

varid = netcdf.defVar(ncid,'dv','NC_DOUBLE',Nc_dimid);
netcdf.putAtt(ncid,varid,'long_name','sea floor depth');
netcdf.putAtt(ncid,varid,'standard_name','sea_floor_depth_below_geoid');
netcdf.putAtt(ncid,varid,'units','m');
netcdf.putAtt(ncid,varid,'positive','down');
netcdf.putAtt(ncid,varid,'coordinates','xv yv');
netcdf.putAtt(ncid,varid,'mesh','suntans_mesh');
netcdf.putAtt(ncid,varid,'location','face');

varid = netcdf.defVar(ncid,'yp','NC_DOUBLE',Np_dimid);
netcdf.putAtt(ncid,varid,'long_name','Northing of 2D mesh node');
netcdf.putAtt(ncid,varid,'standard_name','Northing');

varid = netcdf.defVar(ncid,'yv','NC_DOUBLE',Nc_dimid);
netcdf.putAtt(ncid,varid,'long_name','Northing of 2D mesh face');
netcdf.putAtt(ncid,varid,'standard_name','Northing');

varid = netcdf.defVar(ncid,'cells','NC_INT',[numsides_dimid Nc_dimid ]);
netcdf.putAtt(ncid,varid,'long_name','Maps every face to its corner nodes');
netcdf.putAtt(ncid,varid,'standard_name','face_node_connectivity');

varid = netcdf.defVar(ncid,'nfaces','NC_INT',[Nc_dimid]);
netcdf.putAtt(ncid,varid,'long_name','Number of cell faces');

varid = netcdf.defVar(ncid,'grad','NC_INT',[two_dimid Ne_dimid ]);
netcdf.putAtt(ncid,varid,'long_name','Maps every edge to the two faces it connects');
netcdf.putAtt(ncid,varid,'standard_name','edge_face_connectivity');

varid = netcdf.defVar(ncid,'time','NC_DOUBLE',time_dimid);
netcdf.putAtt(ncid,varid,'units','seconds since 1990-01-01 00:00:00');
netcdf.putAtt(ncid,varid,'long_name','time');
netcdf.defVarFill(ncid,varid,false,999999);

varid = netcdf.defVar(ncid,'eta','NC_DOUBLE',[Nc_dimid time_dimid ]);
netcdf.putAtt(ncid,varid,'units','meters');
netcdf.putAtt(ncid,varid,'long_name','sea surface elevation');
netcdf.putAtt(ncid,varid,'coordinates','time xv yv');
netcdf.defVarFill(ncid,varid,false,999999);

varid = netcdf.defVar(ncid,'uf','NC_DOUBLE',[Ne_dimid Nk_dimid time_dimid]);
netcdf.putAtt(ncid,varid,'units','meters second-1');
netcdf.putAtt(ncid,varid,'long_name','Easward water velocity component');
netcdf.putAtt(ncid,varid,'coordinates','time z_r xv yv');
netcdf.defVarFill(ncid,varid,false,999999);

% varid = netcdf.defVar(ncid,'uc','NC_DOUBLE',[Nc_dimid Nk_dimid time_dimid]);
% netcdf.putAtt(ncid,varid,'units','meters second-1');
% netcdf.putAtt(ncid,varid,'long_name','Easward water velocity component');
% netcdf.putAtt(ncid,varid,'coordinates','time z_r xv yv');
% netcdf.defVarFill(ncid,varid,false,999999);
% 
% varid = netcdf.defVar(ncid,'vc','NC_DOUBLE',[Nc_dimid Nk_dimid time_dimid]);
% netcdf.putAtt(ncid,varid,'units','meters second-1');
% netcdf.putAtt(ncid,varid,'long_name','Northward water velocity component');
% netcdf.putAtt(ncid,varid,'coordinates','time z_r xv yv');
% netcdf.defVarFill(ncid,varid,false,999999);

varid = netcdf.defVar(ncid,'salt','NC_DOUBLE',[Nc_dimid Nk_dimid time_dimid]);
netcdf.putAtt(ncid,varid,'units','ppt');
netcdf.putAtt(ncid,varid,'long_name','Salinity');
netcdf.putAtt(ncid,varid,'coordinates','time z_r xv yv');
netcdf.defVarFill(ncid,varid,false,999999);

varid = netcdf.defVar(ncid,'temp','NC_DOUBLE',[Nc_dimid Nk_dimid time_dimid]);
netcdf.putAtt(ncid,varid,'units','degrees C');
netcdf.putAtt(ncid,varid,'long_name','Water temperature');
netcdf.putAtt(ncid,varid,'coordinates','time z_r xv yv');
netcdf.defVarFill(ncid,varid,false,999999);

varid = netcdf.defVar(ncid,'agec','NC_DOUBLE',[Nc_dimid Nk_dimid time_dimid]);
netcdf.putAtt(ncid,varid,'units','');
netcdf.putAtt(ncid,varid,'long_name','Age concentration');
netcdf.putAtt(ncid,varid,'coordinates','time z_r xv yv');
netcdf.defVarFill(ncid,varid,false,999999);

varid = netcdf.defVar(ncid,'agealpha','NC_DOUBLE',[Nc_dimid Nk_dimid time_dimid]);
netcdf.putAtt(ncid,varid,'units','seconds');
netcdf.putAtt(ncid,varid,'long_name','Age alpha parameter');
netcdf.putAtt(ncid,varid,'coordinates','time z_r xv yv');
netcdf.defVarFill(ncid,varid,false,999999);

varid = netcdf.defVar(ncid,'agesource','NC_DOUBLE',[Nc_dimid Nk_dimid time_dimid]);
netcdf.putAtt(ncid,varid,'units','');
netcdf.putAtt(ncid,varid,'long_name','Age source grid cell (>0 = source)');
netcdf.putAtt(ncid,varid,'coordinates','time z_r xv yv');
netcdf.defVarFill(ncid,varid,false,999999);

% end define variables
netcdf.endDef(ncid);

% write grid variables
varid = netcdf.inqVarID(ncid,'Nk');
netcdf.putVar(ncid,varid,Nk);

varid = netcdf.inqVarID(ncid,'neigh');
netcdf.putVar(ncid,varid,neigh);

varid = netcdf.inqVarID(ncid,'xp');
netcdf.putVar(ncid,varid,xp);

varid = netcdf.inqVarID(ncid,'xv');
netcdf.putVar(ncid,varid,xv);

varid = netcdf.inqVarID(ncid,'z_r');
netcdf.putVar(ncid,varid,z_r);

varid = netcdf.inqVarID(ncid,'latv');
netcdf.putVar(ncid,varid,latv);

varid = netcdf.inqVarID(ncid,'mark');
netcdf.putVar(ncid,varid,mark);

varid = netcdf.inqVarID(ncid,'lonv');
netcdf.putVar(ncid,varid,lonv);

varid = netcdf.inqVarID(ncid,'edges');
netcdf.putVar(ncid,varid,edges);

varid = netcdf.inqVarID(ncid,'dz');
netcdf.putVar(ncid,varid,dz);

varid = netcdf.inqVarID(ncid,'dv');
netcdf.putVar(ncid,varid,dv);

varid = netcdf.inqVarID(ncid,'yp');
netcdf.putVar(ncid,varid,yp);

varid = netcdf.inqVarID(ncid,'yv');
netcdf.putVar(ncid,varid,yv);

varid = netcdf.inqVarID(ncid,'cells');
netcdf.putVar(ncid,varid,cells);

varid = netcdf.inqVarID(ncid,'nfaces');
netcdf.putVar(ncid,varid,nfaces);

varid = netcdf.inqVarID(ncid,'grad');
netcdf.putVar(ncid,varid,grad);



% close file
netcdf.close(ncid);

%% compute fields
disp('computing IC variables')

if initial_ideal %idealized initial conditions
    time = Toffset*86400; % seconds since basetime
        
    % free surface
    eta = zeros(Nc,1);
    
    % velocity
    uc = zeros(Nc,Nkmax,1);
    vc = uc;
    
    % salt
    salt = uc;
    
    % temp
    a1 = 23.3637651923536;
    a2 = 3.18189680534772;
    a3 = -44.1196509264888;
    a4 = 293.118606564064;
    for i=1:Nc
        temp(i,:,1) = a1*exp(-(-z_r+a3)/a4)+a2;
    end
    
    % age parameters
    agec = uc;
    agealpha = uc;
    agesource = uc;    
    
else % realistic initial conditions from NCOM
    
    ncid = netcdf.open(NCOM_file,'NC_NOWRITE');
    fill_value = -99;
    
    varid = netcdf.inqDimID(ncid,'time');
    [~, tcount] = netcdf.inqDim(ncid,varid);
    tstart=0; % starts at 0
    
    % time
    varname = netcdf.inqVar(ncid,0);
    varid = netcdf.inqVarID(ncid,varname);
    temp_time_all = netcdf.getVar(ncid,varid,tstart,tcount,'double');
    NCOM_basetime = datenum(netcdf.getAtt(ncid,varid,'start'));
    time_step = netcdf.getAtt(ncid,varid,'units');
    if strcmp(time_step,'hour')
        time_step = 24;
    end
    temp_mtime_all = temp_time_all/time_step+NCOM_basetime;

    % lon
    varname = netcdf.inqVar(ncid,1);
    varid = netcdf.inqVarID(ncid,varname);
    temp_lon = netcdf.getVar(ncid,varid,'double');

    % lat
    varname = netcdf.inqVar(ncid,2);
    varid = netcdf.inqVarID(ncid,varname);
    temp_lat = netcdf.getVar(ncid,varid,'double');

    %depth
    varname = netcdf.inqVar(ncid,3);
    varid = netcdf.inqVarID(ncid,varname);
    temp_depth = netcdf.getVar(ncid,varid,'double');
    temp_z = -temp_depth;

    temp_xcount = length(temp_lon);
    temp_ycount = length(temp_lat);
    temp_zcount = length(temp_depth);
    
    [temp_LAT,temp_LON]=meshgrid(temp_lat,temp_lon);    
    [temp_Y,temp_X,Zone,lcm]=ell2utm(deg2rad(temp_LAT),deg2rad(temp_LON));
    lcm = nanmedian(reshape(lcm(Zone==GRID.Zone),[],1));
    [temp_Y,temp_X,~,~]=ell2utm(deg2rad(temp_LAT),deg2rad(temp_LON),lcm);
    
    temp_dy = max(max(abs(diff(temp_Y,[],2))));
    temp_dx = max(max(abs(diff(temp_X,[],1))));
    
    % get buffer points, ccw from ll, dx from model    
    x1 = min(xv)-10*temp_dx;
    x2 = max(xv)+10*temp_dx;    
    y1 = min(yv)-10*temp_dy;
    y2 = max(yv)+10*temp_dy;    
    xbuffer = [x1 x2 x2 x1 x1];
    ybuffer = [y1 y1 y2 y2 y1];
    
    % find spatial indices (NC files start at 0)
    indx = inpolygon(temp_X,temp_Y,xbuffer,ybuffer);
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
    zcount = length(temp_depth);
    
    % find time index closest to starttime
    [~,tindx] = min(abs(mtime_start-temp_mtime_all));
    tindx = tindx-1; %netcdf indexing
    %% 
    % reload time, lat, lon data only in area of interest
    
    % time
    varname = netcdf.inqVar(ncid,0);
    varid = netcdf.inqVarID(ncid,varname);
    temp_time = netcdf.getVar(ncid,varid,tindx,1,'double');
    temp_mtime =  temp_time/24+NCOM_basetime;
    time = (temp_mtime-mtime_base)*86400; % time in s since basetime
  
    % lon
    varname = netcdf.inqVar(ncid,1);
    varid = netcdf.inqVarID(ncid,varname);
    temp_lon = netcdf.getVar(ncid,varid,[xstart],[xcount],'double');

    % lat
    varname = netcdf.inqVar(ncid,2);
    varid = netcdf.inqVarID(ncid,varname);
    temp_lat = netcdf.getVar(ncid,varid,ystart,ycount,'double');
    
    [temp_LAT,temp_LON]=meshgrid(temp_lat,temp_lon);    
    [temp_Y,temp_X,Zone,lcm]=ell2utm(deg2rad(temp_LAT),deg2rad(temp_LON),lcm);
       
    % temperature
    disp(['interpolating IC from NCOM netcdf ' datestr(temp_mtime)])
    fprintf('eta...')
    varid = netcdf.inqVarID(ncid,'ssh');
    temp_ssh = netcdf.getVar(ncid,varid,[xstart ystart tindx],[xcount ycount 1],'double');
    temp_ssh(temp_ssh==fill_value)=nan;
    ydata = reshape(temp_ssh,[],1);
    xdata = reshape(temp_X,[],1);
    xdata2 = reshape(temp_Y,[],1);
    indxe = ~isnan(ydata); % in case any dry cells
    [px]=polyfit(xdata(indxe),ydata(indxe),1);
    [py]=polyfit(xdata2(indxe),ydata(indxe),1);
    etatide = px(1)*temp_X + py(1)*temp_Y + px(2);
    etawave = temp_ssh - etatide;
    etawave = etawave - nanmean(nanmean(etawave));
    % remove the wave component for stability
    ssh_smooth = temp_ssh-etawave;     
    % interp to grid
    % use original ssh data
    [eta] = interpNCOM_parallel(temp_X,temp_Y,temp_z,temp_ssh,xv,yv,z_r);
    % use smoothed ssh
%     [eta] = interpNCOM_parallel(temp_X,temp_Y,temp_z,ssh_smooth,xv,yv,z_r);
    % set I/C to flat,
    if etainterp==0
       eta = 0*eta; 
    end
    
    fprintf('T...')
    varid = netcdf.inqVarID(ncid,'temperature');
    NC_load = netcdf.getVar(ncid,varid,[xstart ystart 0 tindx],[xcount ycount temp_zcount 1],'double');
    NC_load(NC_load==fill_value)=nan;
    % interp to grid
    [temp] = interpNCOM_parallel(temp_X,temp_Y,temp_z,NC_load,xv,yv,z_r);
     
     % salinity
    fprintf('S...')
    varid = netcdf.inqVarID(ncid,'salinity');
    NC_load = netcdf.getVar(ncid,varid,[xstart ystart 0 tindx],[xcount ycount temp_zcount 1],'double');
    NC_load(NC_load==fill_value)=nan;   
    % interp to grid
    [salt] = interpNCOM_parallel(temp_X,temp_Y,temp_z,NC_load,xv,yv,z_r);
        
    % u
    fprintf('u...')
    varid = netcdf.inqVarID(ncid,'u-velocity');
    NC_load = netcdf.getVar(ncid,varid,[xstart ystart 0 tindx],[xcount ycount temp_zcount 1],'double');
    NC_load(NC_load==fill_value)=nan;   
    % interp to grid
    [ue] = interpNCOM_parallel(temp_X,temp_Y,temp_z,NC_load,xe,ye,z_r);
    
    % v
    fprintf('v...')
    varid = netcdf.inqVarID(ncid,'v-velocity');
    NC_load = netcdf.getVar(ncid,varid,[xstart ystart 0 tindx],[xcount ycount temp_zcount 1],'double');
    NC_load(NC_load==fill_value)=nan;   
    % interp to grid
    [ve] = interpNCOM_parallel(temp_X,temp_Y,temp_z,NC_load,xe,ye,z_r);
    
    % compute u at cell faces
    fprintf('uface...')
    uf = N1.*ue + N2.*ve;
    
    % extrapolate nans in case of mismatch in values
%     eta = naninterp(eta);
%     temp = naninterp(temp);
%     salt = naninterp(salt);
    
    % set u,v,eta to zero for IC
%     fprintf('set u=v=0...')
%     eta = zeros(Nc,1);
    uc = 0*temp;
    vc = 0*temp;    
    agec = 0*temp;
    agealpha = 0*temp;
    agesource = 0*temp;  
    fprintf('done!\n')
    
    
    netcdf.close(ncid);
    
end





%% write variables
disp('writing variables to netcdf IC file')

ncid = netcdf.open([datadir init_filename],'WRITE');

varid = netcdf.inqVarID(ncid,'time');
netcdf.putVar(ncid, varid,0,1, time);

varid = netcdf.inqVarID(ncid,'eta');
netcdf.putVar(ncid, varid,eta);

varid = netcdf.inqVarID(ncid,'uf');
netcdf.putVar(ncid, varid,uf);

% varid = netcdf.inqVarID(ncid,'uc');
% netcdf.putVar(ncid, varid,uc);
% 
% varid = netcdf.inqVarID(ncid,'vc');
% netcdf.putVar(ncid, varid,vc);

varid = netcdf.inqVarID(ncid,'salt');
netcdf.putVar(ncid, varid,salt);

varid = netcdf.inqVarID(ncid,'temp');
netcdf.putVar(ncid, varid,temp);

varid = netcdf.inqVarID(ncid,'agec');
netcdf.putVar(ncid, varid,agec);

varid = netcdf.inqVarID(ncid,'agealpha');
netcdf.putVar(ncid, varid,agealpha);

varid = netcdf.inqVarID(ncid,'agesource');
netcdf.putVar(ncid, varid,agesource);

netcdf.close(ncid);


%% plot IC conditions

xplot = reshape(xv,GRID.Nx,GRID.Ny)/1000;
yplot = reshape(yv,GRID.Nx,GRID.Ny)/1000;
xeplot = xe/1000;
yeplot = ye/1000;
etaplot = reshape(eta,GRID.Nx,GRID.Ny);
Tplot = reshape(temp(:,1),GRID.Nx,GRID.Ny);
Splot = reshape(salt(:,1),GRID.Nx,GRID.Ny);
UplotM = squeeze(sqrt(uf(:,1).^2));
% Uplot = reshape(uc(:,1),GRID.Nx,GRID.Ny);
% Vplot = reshape(vc(:,1),GRID.Nx,GRID.Ny);
% UplotM = sqrt(Uplot.^2+Vplot.^2);

figure(3452423)
subplot(2,2,1)
if GRID.Ny>2
pcolorjw(xplot,yplot,etaplot)
cb=colorbar;
ylabel(cb,'$\eta$ (m)','interpreter','latex')
ylabel('y (km)')
else
    plot(xplot,etaplot,'.')
    ylabel('$\eta$ (m)','interpreter','latex')
end

if exist('NC')
    title(['IC netcdf ' datestr(temp_mtime)])
else
    title(['IC ideal ' datestr(mtime_start)])
end
    

subplot(2,2,2)
if GRID.Ny>2
scatter(xeplot,yeplot,5,UplotM)
cb=colorbar;
ylabel(cb,'$|u|_{surf}$ (m/s)','interpreter','latex')
else
    plot(xeplot,UplotM,'.')
    ylabel('$|u|_{surf}$ (m/s)','interpreter','latex')
end

subplot(2,2,3)
if GRID.Ny>2
pcolorjw(xplot,yplot,Tplot)
cb=colorbar;
ylabel(cb,'$T_{surf}$ (C)','interpreter','latex')
ylabel('y (km)')
else
    plot(xplot,Tplot)
    ylabel('$T_{surf}$ (C)','interpreter','latex');
end

subplot(2,2,4)
if GRID.Ny>2
pcolorjw(xplot,yplot,Splot)
cb=colorbar;
ylabel(cb,'$S_{surf}$ (psu)','interpreter','latex')
else
    plot(xplot,Splot)
    ylabel('$S_{surf}$ (psu)','interpreter','latex')
end
xlabel('x (km)')

print -djpeg -r300 figure_IC_NCOM
close


delete(poolobj)






