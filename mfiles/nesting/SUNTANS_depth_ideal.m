function [] = SUNTANS_depth_ideal(datadir,D0,shape,bathy_dir,x0,y0,theta,dtrend,dmin,extreme)

% datadir='../rundata';
% H=3000;
if nargin<9
    dmin=0;
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

points = load([datadir,'/points.dat']);
Np = size(points,1);
xp = points(:,1);
yp = points(:,2);

ed = load([datadir,'/edges.dat']);
Ne = size(ed,1);
edges = ed(:,1:2);
mark = ed(:,3);
grad =ed(:,4:5);

load('SUNTANS_grid.mat','Nx','L','W','Ny','Zone')
dx = L/Nx;
dy = W/Ny;

sponge_distance =  getvalue([datadir,'/suntans.dat'],'sponge_distance');

%% create z at cell points xv yv
if strcmp(shape,'tanh')

    Ls = 30000/198000*(max(xv)-min(xv));
    xmid = nanmean(xv);
    Ds = D0*200/3000;

    Depth = D0 - 0.5*(D0-Ds)*(1+tanh(4*(xv-xmid)/Ls));

elseif strcmp(shape,'const')
    
    Depth = D0+0*xv;
    
elseif strcmp(shape,'witch')
    xmid = nanmean(xv);
    h0 = D0/10;
    b = 4*D0;
    
    Depth = D0 - h0.*(1+(xv-xmid).^2/b^2).^(-1);
    
elseif strcmp(shape,'sine1D')
    load('./SUNTANS_grid.mat')
    a = 0.9*1.188;
    Depth = D0 + a*cos(2*pi*2/L*xv);
    
elseif strcmp(shape,'sine2D')
    load('./SUNTANS_grid.mat')
    a = 1.0;
    Depth = D0 + a*cos(2*pi*2/L*xv).*cos(2*pi*1/W*yv);
    
elseif strcmp(shape,'square')
    load('./SUNTANS_grid.mat')
    Depth = D0 + 0.*xv;
    xmid = nanmean(xv);
    ymid = nanmean(yv);
    
    hs = D0/2;
    Ws = W/3;
    
    indx = xv >= xmid-Ws/2 & xv <= xmid+Ws/2 &...
        yv >= ymid-Ws/2 & yv <= ymid+Ws/2;
    
    Depth (indx) = hs;
    
elseif strcmp(shape,'one_cell')
    load('./SUNTANS_grid.mat')
    Depth = D0 + 0.*xv;
    xmid = nanmean(xv);
    ymid = nanmean(yv);
    Nk=4;
    dz = D0/Nk;
    dx=L/Nx;
    dy=W/Ny;
    
    hs = 1*dz;
    
    % this is for 4 grid cells
%     indx = xv > xmid-dx & xv <= xmid+dx &...
%         yv > ymid-dy & yv <= ymid+dy;
    
    % this is for one grid cell
    r = sqrt((xmid-xv).^2 + (ymid-yv).^2);
    [~,indx] = min(r);
   
    Depth (indx) = D0-hs;    
    
elseif strcmp(shape,'angle')
    load('./SUNTANS_grid.mat')
    Depth = D0 + 0.*xv;
    xmid = nanmean(xv);
    ymid = nanmean(yv);
    
    hs = D0/2;
    Ws = W/2;
    
    % left side
    indx = xv >= xmid-hs & xv <= xmid &...
        yv >= ymid-Ws/2 & yv <= ymid+Ws/2;
    Depth (indx) = hs-(xv(indx)-xmid);
        
    % right side
    indx = xv > xmid & xv <= xmid+hs &...
        yv >= ymid-Ws/2 & yv <= ymid+Ws/2;
    Depth (indx) = hs+(xv(indx)-xmid);
    
elseif strcmp(shape,'cylinder')
    load('./SUNTANS_grid.mat')
    Depth = D0 + 0.*xv;
    xmid = nanmean(xv);
    ymid = nanmean(yv);
    
    hs = D0/2;
    Ws = W/3;
    
    r = sqrt((xv-xmid).^2+(yv-ymid).^2);
    
    indx = r<=Ws/2;
    
    Depth (indx) = hs;
    
elseif strcmp(shape,'real')
    method='cubic';
    
    % load bathy data
    DATA = load(bathy_dir);
    x = reshape(DATA.X,[],1);
    y = reshape(DATA.Y,[],1);
    z = reshape(DATA.Z,[],1);
    
    x_utm = reshape(DATA.X_UTM,[],1);
    y_utm = reshape(DATA.Y_UTM,[],1);
    % put xv, yv in reference frame of bathym
    
    % rotate xv and yv to absolute ref frame
    xv_ref = xv*cosd(theta)-yv*sind(theta);
    yv_ref = xv*sind(theta)+yv*cosd(theta);
    % locate lower left corner at x0,y0
    xv_ref = xv_ref + x0;
    yv_ref = yv_ref + y0;
    % compute coords in UTM for referencing
    xv_UTM = xv_ref + (x_utm(1)-x(1));
    yv_UTM = yv_ref + (y_utm(1)-y(1));
    
    L = max(xv)-min(xv);
    W = max(yv)-min(yv);
    
    % find points nearby,
    indx = x < max(xv_ref)+L/2 & x > min(xv_ref)-L/2 &...
        y < max(yv_ref)+W/2 & y > min(yv_ref)-W/2 &...
        ~isnan(z);  
    
    
    
%     F = scatteredInterpolant(x(indx),y(indx),z(indx));
%     zv = F(xv_ref,yv_ref);
    
    zv = griddata(x(indx),y(indx),z(indx),xv_ref,yv_ref,method);
    
    zv = zv - nanmean(zv);
     
    
    % remove high and low extremes
    if nargin<10
        extreme=0;
    end
    if extreme        
        indx = zv> quantile(zv,0.02) & zv < quantile(zv,0.98);
        F = scatteredInterpolant(xv(indx),yv(indx),zv(indx));
        zv(~indx) = F(xv(~indx),yv(~indx));
    end
        
    if dtrend
        % detrend in x
        p=polyfit(xv,zv,1);
        zv_fit = polyval(p,xv);
        zv = zv - zv_fit;
        
        % detrend in y
        p=polyfit(yv,zv,1);
        zv_fit = polyval(p,yv);
        zv = zv - zv_fit;
    end
    %
    
    alpha=2;
    fx = 0.5*(tanh(xv-alpha-min(min(xv)))+...
    tanh(-xv-alpha+max(max(xv))));

    Zgrid = zv.*fx;
    Zgrid = Zgrid - nanmean(Zgrid);
    
    % add in depth
    Depth = -Zgrid+D0;
    
    % check for D<Dmin
    Depth(Depth<dmin)=dmin;
    
    % find corners
    x1 = min(xv);
    x2 = max(xv);
    y1 = min(yv);
    y2 = max(yv);
    
    xcorn = [x1 x2 x2 x1 x1];
    ycorn = [y1 y1 y2 y2 y1];
    
    % rotate xcorn ycorn to absolute ref frame
    xcorn_ref = xcorn*cosd(theta)-ycorn*sind(theta);
    ycorn_ref = xcorn*sind(theta)+ycorn*cosd(theta);
    % locate lower left corner at x0,y0
    xcorn_ref = xcorn_ref + x0;
    ycorn_ref = ycorn_ref + y0;
    
    
    % compute stats on bathym
    BATH.x = reshape(xv,Nx,Ny);
    BATH.y = reshape(yv,Nx,Ny);
    
    BATH.x_UTM = reshape(xv_UTM,Nx,Ny);
    BATH.y_UTM = reshape(yv_UTM,Nx,Ny);
    
    BATH.z = z;
    BATH.Z = reshape(-Depth,Nx,Ny);
    BATH.D0 = D0;

    BATH.dx = BATH.x(2,2)-BATH.x(1,1);
    BATH.dy = BATH.y(2,2)-BATH.y(1,1);

    BATH.dzdx = nan+BATH.x;
    BATH.dzdx(2:end,:) = (BATH.Z(2:end,:) - BATH.Z(1:end-1,:))/BATH.dx;

    BATH.dzdy = nan+BATH.y;
    BATH.dzdy(:,2:end) = (BATH.Z(:,2:end) - BATH.Z(:,1:end-1))/BATH.dy;

%     BATH.Ub = DRAG.Ub;
%     BATH.nu = 1e-6;

    BATH.hrms = sqrt(8*nanmean(reshape((BATH.Z+BATH.D0).^2,[],1)));

    BATH.dzdx_rms = sqrt(nanmean(reshape((BATH.dzdx).^2,[],1)));
    BATH.dzdy_rms = sqrt(nanmean(reshape((BATH.dzdy).^2,[],1)));


    % fourier transform of depth and slope
    for i=1:size(BATH.x,2) % go through each y coord
        % (space, freq)
       [~,BATH.Fzx(i,:),BATH.Sx]=low_pass_filter(BATH.Z(:,i)+BATH.D0,BATH.dx,BATH.dx);
       [~,BATH.Fdzdx(i,:),~]=low_pass_filter(BATH.dzdx(:,i),BATH.dx,BATH.dx);
    end
    if size(BATH.x,2)>3
    for i=1:size(BATH.x,1) % go through each x coord
        % (space, freq)
       [~,BATH.Fzy(i,:),BATH.Sy]=low_pass_filter(BATH.Z(i,:)+BATH.D0,BATH.dy,BATH.dy);
       [~,BATH.Fdzdy(i,:),~]=low_pass_filter(BATH.dzdy(i,:),BATH.dy,BATH.dy);
    end
    end
    
    BATH.dfx = 1./abs(BATH.x(2,2)-BATH.x(1,1));
    BATH.dfy = 1./abs(BATH.y(2,2)-BATH.y(1,1));
    
    % get energy spectral density
    BATH.Szx = 2*BATH.Fzx.*conj(BATH.Fzx);
    BATH.Sdzdx = 2*BATH.Fdzdx.*conj(BATH.Fdzdx);
    
    if size(BATH.x,2)>3
    BATH.Szy = 2*BATH.Fzy.*conj(BATH.Fzy);
    BATH.Sdzdy = 2*BATH.Fdzdy.*conj(BATH.Fdzdy);
    end

    % take mean in space x
    BATH.Szx_bar = nanmean(BATH.Szx,1);
    BATH.Sdzdx_bar = nanmean(BATH.Sdzdx,1);
    % take mean in space y
    if size(BATH.x,2)>3
        BATH.Szy_bar = nanmean(BATH.Szy,1);
        BATH.Sdzdy_bar = nanmean(BATH.Sdzdy,1);
    else
       BATH.Szy_bar=nan;
       BATH.Sdzdy_bar = nan;
    end
    
    % find wavelength of centroid
    BATH.Lx = 1./BATH.Sx;    
    BATH.Lxmax_avg = sum(BATH.Sdzdx_bar)./sum(BATH.Sx.*BATH.Sdzdx_bar);
%     BATH.Lxmax_peak = sqrt(sum(BATH.Sx.^2.*BATH.Sdzdx_bar)./sum(BATH.Sx.^4.*BATH.Sdzdx_bar));
    [~,indx]=max(abs(BATH.Sdzdx_bar));
    BATH.Lxmax_peak = BATH.Lx(indx);

    if size(BATH.x,2)>3
        BATH.Ly = 1./BATH.Sy;    
        BATH.Lymax_avg = sum(BATH.Sdzdy_bar)./sum(BATH.Sy.*BATH.Sdzdy_bar);
    %     BATH.Lxmax_peak = sqrt(sum(BATH.Sx.^2.*BATH.Sdzdx_bar)./sum(BATH.Sx.^4.*BATH.Sdzdx_bar));
        [~,indx]=max(abs(BATH.Sdzdy_bar));
        BATH.Lymax_peak = BATH.Ly(indx);
    else
        BATH.Ly=nan;
        BATH.Lymax_avg=nan;
        BATH.Lymax_peak = nan;
    end

    % compute stats
    BATH.h_H = BATH.hrms/BATH.D0;
%     BATH.Re_h = BATH.Ub*BATH.hrms/BATH.nu;
    BATH.h_lambdax = sqrt(2)/pi*BATH.dzdx_rms;
    if size(BATH.x,2)>3
        BATH.h_lambday = sqrt(2)/pi*BATH.dzdy_rms;
    else
        BATH.h_lambday=0;
    end
    
    save('SUNTANS_grid.mat','-append','BATH','bathy_dir','theta');
    
elseif strcmp(shape,'real_interp_only')
    Dmax=D0;
    
    % load bathy data
    DATA = load(bathy_dir);
    x = reshape(DATA.X,[],1);
    y = reshape(DATA.Y,[],1);
    z = reshape(DATA.Z,[],1);
    
    x_utm = reshape(DATA.X_UTM,[],1);
    y_utm = reshape(DATA.Y_UTM,[],1);
    % put xv, yv in reference frame of bathym
    
    % rotate xv and yv to absolute ref frame
    xv_ref = xv*cosd(theta)-yv*sind(theta);
    yv_ref = xv*sind(theta)+yv*cosd(theta);
    % locate lower left corner at x0,y0
    xv_ref = xv_ref + x0;
    yv_ref = yv_ref + y0;
    % compute coords in UTM for referencing
    xv_UTM = xv_ref + (x_utm(1)-x(1));
    yv_UTM = yv_ref + (y_utm(1)-y(1));
    
    L = max(xv)-min(xv);
    W = max(yv)-min(yv);
    
    
    F = scatteredInterpolant(x,y,z);
    
    zv = F(xv_ref,yv_ref);

    % add in depth
    Depth = -zv;
    Depth(Depth>Dmax)=Dmax;
    Depth(Depth<dmin)=dmin;
    

    % find corners
    x1 = min(xv);
    x2 = max(xv);
    y1 = min(yv);
    y2 = max(yv);
    
    xcorn = [x1 x2 x2 x1 x1];
    ycorn = [y1 y1 y2 y2 y1];
    
    % rotate xcorn ycorn to absolute ref frame
    xcorn_ref = xcorn*cosd(theta)-ycorn*sind(theta);
    ycorn_ref = xcorn*sind(theta)+ycorn*cosd(theta);
    % locate lower left corner at x0,y0
    xcorn_ref = xcorn_ref + x0;
    ycorn_ref = ycorn_ref + y0;
    
    
    % compute stats on bathym
    BATH.x = reshape(xv,Nx,Ny);
    BATH.y = reshape(yv,Nx,Ny);
    
    BATH.x_UTM = reshape(xv_UTM,Nx,Ny);
    BATH.y_UTM = reshape(yv_UTM,Nx,Ny);
    
    BATH.z = z;
    BATH.Z = reshape(-Depth,Nx,Ny);
    BATH.D0 = D0;

    BATH.dx = dx;
    BATH.dy = dy;

   
    save('SUNTANS_grid.mat','-append','BATH','bathy_dir','theta');
    
elseif strcmp(shape,'NCOM')
    
    file = dir(bathy_dir);
    ncfile=regexp(file.name, regexptranslate('wildcard','*.nc'));
    matfile=regexp(file.name, regexptranslate('wildcard','*.mat'));
    if ncfile % bathymetry is from NCOM netcdf file
        
    % load in NC data
%     bathy_dir = ['C:\Users\jsrogers\Desktop\Transfer\Modeling\SUNTANS\NCOM_SCS\lzsnfs-201312.nc'];
    ncid = netcdf.open(bathy_dir,'NC_NOWRITE');
    fill_value = -99;
    
    varid = netcdf.inqDimID(ncid,'time');
    [~, tcount] = netcdf.inqDim(ncid,varid);
    tstart=0; % starts at 0
    
    % time
    varname = netcdf.inqVar(ncid,0);
    varid = netcdf.inqVarID(ncid,varname);
    NC.time = netcdf.getVar(ncid,varid,tstart,tcount,'double');

    % long
    varname = netcdf.inqVar(ncid,1);
    varid = netcdf.inqVarID(ncid,varname);
    NC.lon = netcdf.getVar(ncid,varid,'double');

    % lat
    varname = netcdf.inqVar(ncid,2);
    varid = netcdf.inqVarID(ncid,varname);
    NC.lat = netcdf.getVar(ncid,varid,'double');

    %depth
    varname = netcdf.inqVar(ncid,3);
    varid = netcdf.inqVarID(ncid,varname);
    NC.depth = netcdf.getVar(ncid,varid,'double');
    
    % get depth on bottom of cell face
    NC.z = -NC.depth;
    NC.zface = [0; 0.5*(NC.z(1:end-1)+NC.z(2:end))];
    NC.zface = [NC.zface; NC.z(end)+(NC.z(end)-NC.zface(end))];
    NC.depth_bot_face = -NC.zface(2:end);
    NC.dz = diff(NC.zface);

    NC.xcount = length(NC.lon);
    NC.ycount = length(NC.lat);
    NC.zcount = length(NC.depth);
    
    % u
    varname = netcdf.inqVar(ncid,7);
    varid = netcdf.inqVarID(ncid,varname);
    NC.u = netcdf.getVar(ncid,varid,[0 0 0 0],[NC.xcount NC.ycount NC.zcount 1],'double');
    NC.u(NC.u==fill_value)=nan;
    
    netcdf.close(ncid);
    
    [NC.LAT,NC.LON]=meshgrid(NC.lat,NC.lon);

    % get full depth grid
    for i=1:length(NC.lon)
        for j=1:length(NC.lat)
            indx = find(~isnan(squeeze(NC.u(i,j,:,1))));
            if ~isempty(indx)
                NC.Depth(i,j) = NC.depth_bot_face(indx(end));
            else
                NC.Depth(i,j)=nan;
            end
        end
    end
    
    [NC.Y,NC.X,Zone,lcm]=ell2utm(deg2rad(NC.LAT),deg2rad(NC.LON));
    Zone = median(reshape(Zone,[],1));
    
%     NC.x = reshape(NC.X,[],1);
%     NC.y = reshape(NC.Y,[],1);
%     NC.z = -reshape(NC.Depth,[],1);
    
    NC.dy = max(max(abs(diff(NC.Y,[],2))));
    NC.dx = max(max(abs(diff(NC.X,[],1))));
    
    % get buffer points, ccw from ll, dx from model    
    x1 = min(xv)-2*NC.dx;
    x2 = max(xv)+2*NC.dx;    
    y1 = min(yv)-2*NC.dy;
    y2 = max(yv)+2*NC.dy;    
    xbuffer = [x1 x2 x2 x1 x1];
    ybuffer = [y1 y1 y2 y2 y1];
    
    % this is the index to use
    indx = inpolygon(NC.X,NC.Y,xbuffer,ybuffer);
    % rows are x, columns are y
    xsum = find(sum(indx,2)>0);
    xstart=xsum(1);
    xend =xsum(end);
    % rows are x, columns are y
    ysum = find(sum(indx,1)>0);
    ystart=ysum(1);
    yend =ysum(end);
    
    
    Dmax=max(max(NC.Depth));
        
    % create surface
    x = reshape(NC.X(indx),[],1);
    y = reshape(NC.Y(indx),[],1);
    z = reshape(-NC.Depth(indx),[],1); 
    
   BATH.NC = NC;
    
    elseif matfile % bathy is from .mat file
        MAT = load(bathy_dir);        
        [MAT.y,MAT.x,MAT.Zone,MAT.lcm]=ell2utm(deg2rad(MAT.lat),deg2rad(MAT.lon));
        lcm = nanmedian(MAT.lcm(MAT.Zone==Zone)); % get lcm for specified UTM zone;
        % now recompute with constant lcm
        [MAT.y,MAT.x,MAT.Zone,MAT.lcm]=ell2utm(deg2rad(MAT.lat),deg2rad(MAT.lon),lcm);
       
        Zone = median(reshape(Zone,[],1));
        % select buffer region
        offset = max([(max(xv) - min(xv))/10, (max(yv) - min(yv))/10]);
        indx = (MAT.x) > (min(xv)-offset) &...
               (MAT.x) < (max(xv)+offset) &...
               (MAT.y) > (min(yv)-offset) &...
               (MAT.y) < (max(yv)+offset);
           
        x = reshape(MAT.x(indx),[],1);
        y = reshape(MAT.y(indx),[],1);
        z = reshape(MAT.z(indx),[],1);      
        
        BATH.MAT = MAT;
        
    
    else
        disp('no valid bathymetry file type entered')
        return
    end
    
    
    F = scatteredInterpolant(x,y,z);
    
    zv = F(xv,yv);

    % add in depth
    Depth = -zv;
    
    Depth(Depth<D0)=D0;
    
    Depth(Depth>dmin)=dmin;
    
    % compute stats on bathym
    BATH.x = reshape(xv,Nx,Ny);
    BATH.y = reshape(yv,Nx,Ny);
       
    BATH.z = z;
    BATH.Z = reshape(-Depth,Nx,Ny);
    BATH.D0 = D0;

    BATH.dx = dx;
    BATH.dy = dy;
    
    % add in sponge distances for easy plotting
    % trasnform back into non-rotated frame
    load('./SUNTANS_grid.mat','x0','y0','theta','xv','yv')
    xtemp = (xv-x0)*cosd(-theta)-(yv-y0)*sind(-theta);
    ytemp = (xv-x0)*sind(-theta)+(yv-y0)*cosd(-theta);
    % get points
    xlims = [min(min(xtemp)) max(max(xtemp))];
    ylims = [min(min(ytemp)) max(max(ytemp))];
    % get sponge points
    x_sponge = sponge_distance;
    xpt = [xlims(1)+x_sponge,xlims(2)-x_sponge,xlims(2)-x_sponge,xlims(1)+x_sponge,xlims(1)+x_sponge];
    ypt = [ylims(1)+x_sponge,ylims(1)+x_sponge,ylims(2)-x_sponge,ylims(2)-x_sponge,ylims(1)+x_sponge];
    % transorm to rotated frame
    xpt2 = xpt*cosd(theta)-ypt*sind(theta);
    ypt2 = xpt*sind(theta)+ypt*cosd(theta);
    BATH.sponge_x = xpt2+x0;
    BATH.sponge_y = ypt2+y0;
    
    % get corner points, LLC CCW
    xpt = [xlims(1) xlims(2) xlims(2) xlims(1) xlims(1)];
    ypt = [ylims(1) ylims(1) ylims(2) ylims(2) ylims(1)];
    xpt2 = xpt*cosd(theta)-ypt*sind(theta);
    ypt2 = xpt*sind(theta)+ypt*cosd(theta);
    BATH.corner_x = xpt2+x0;
    BATH.corner_y = ypt2+y0;
    
    % get midpoints of edges, CCW from S
    xpt = [mean(xlims) max(xlims) mean(xlims) min(xlims)];
    ypt = [min(ylims) mean(ylims) max(ylims) mean(ylims)];
    xpt2 = xpt*cosd(theta)-ypt*sind(theta);
    ypt2 = xpt*sind(theta)+ypt*cosd(theta);
    BATH.mid_x = xpt2+x0;
    BATH.mid_y = ypt2+y0;
    
    save('SUNTANS_grid.mat','BATH','bathy_dir','theta','Zone','-append');

end


%% save to file

depthoutput = [xv yv Depth];

depth_file = [datadir,'/depth.dat'];
depthf = fopen(depth_file,'w');
fprintf(depthf, '%12.10e %12.10e %22.20e\n', depthoutput');
status = fclose(depthf);  
%%
save('SUNTANS_grid.mat','Depth','xv','yv','-append');

%% save grid file

plot_grid(datadir)

%% plot depth
close all

subplot(2,1,1)
plot(xv,-Depth,'.')
xlabel('x (m)')
ylabel('z (m)')

subplot(2,1,2)
scatter(xv,yv,5,-Depth)
cb=colorbar;
colormap jet
ylabel(cb,'z (m)')
ylabel('y (m)')
xlabel('x (m)')

print -djpeg -r300 figure_depth
close

%%

if exist('DATA','var')
    
    subplot(1,2,1)
   
    plot(xcorn_ref,ycorn_ref,'-k')
    hold on
    
       
    pcolorjw(DATA.X,DATA.Y,DATA.Z)
    colorbar
    colormap jet
    
    hold on
    plot(xcorn_ref,ycorn_ref,'-k')
    axis equal
    axis manual
    if isfield(DATA,'SITE')
        hold on
        plot(DATA.SITE.x,DATA.SITE.y,'xk')
        hold on
        text(DATA.SITE.x,DATA.SITE.y,DATA.SITE.name)
        
    end
    
    ylabel('y')
    xlabel('x')
    
    
    subplot(1,2,2)
   
    plot(xcorn_ref,ycorn_ref,'-k')
    hold on
    axis equal
    axis manual
    
    pcolorjw(DATA.X,DATA.Y,DATA.Z)
    cb=colorbar;
    caxis([min(-Depth) max(-Depth)])
    colormap jet
    
    hold on
    plot(xcorn_ref,ycorn_ref,'-k')
    
        if isfield(DATA,'SITE')
        hold on
        plot(DATA.SITE.x,DATA.SITE.y,'xk')
        hold on
        text(DATA.SITE.x,DATA.SITE.y,DATA.SITE.name)
        
    end
    
    ylabel('y')
    xlabel('x')
    ylabel(cb,'D')
    
    print -djpeg -r300 figure_bathy
    
    
    %% spectral space of bottom
    if isfield(BATH,'Sx')
    
    clear lgd 

    figure(451635)

    subplot(1,2,1)    
    x = 1./BATH.Sx;
    y = BATH.Szx_bar;
    if size(BATH.x,2)>3
        xx = 1./BATH.Sy;
        yy = BATH.Szy_bar;
    else
        xx = x;
        yy = nan+xx;
    end
    loglog(x,y,xx,yy)

    text(0.1,0.9,['h_{rms}= ' num2str(BATH.hrms) ' m'],'units','normalized');


    subplot(1,2,2)
    x = 1./BATH.Sx;
    y = BATH.Sdzdx_bar;
    if size(BATH.x,2)>3
        xx = 1./BATH.Sy;
        yy = BATH.Sdzdy_bar;
    else
        xx=x;
        yy = nan+xx;
    end
    loglog(x,y,xx,yy)
    hold on
    % 
    text(0.1,0.9,['dz/d(x,y)_{rms}= (' num2str(BATH.dzdx_rms) ',' num2str(BATH.dzdy_rms) ') m'],'units','normalized');
    text(0.1,0.8,['\lambda_{bar}(x,y)= (' num2str(BATH.Lxmax_avg) ',' num2str(BATH.Lymax_avg) ') m'],'units','normalized');

    subplot(1,2,1)
    ylabel('S_{z_{b}}')
    xlabel('\lambda (m)')
    % ylim([1e-7 1e-2])
    text(0.02, 1.05,'(a) depth','units','normalized')

    subplot(1,2,2)
    ylabel('S_{dz_b/dx}')
    xlabel('\lambda (m)')
    % ylim([1e-4 1e-2])
    text(0.02, 1.05,'(b) slope','units','normalized')

    legend('x','y','location','se')

    % get(gcf,'position')
    set(gcf,'position',[680   712   455   266])


    print -djpeg -r300 figure_bathy_spectra_bottom
    end
end
    %%
    if exist('BATH','var')
        load('./SUNTANS_grid_quadgrid.mat','xe','ye','mark')
        figure(23423)
        if isfield(BATH,'NC')
            
        
        pcolorjw(BATH.NC.X/1000,BATH.NC.Y/1000,-BATH.NC.Depth)
        cb=colorbar;
        ylabel(cb,'z (m)')
        
        hold on
        indx = mark==3 | mark ==2;
        plot(xe(indx),ye(indx),'k.')
%         x1 = min(xv);
%         x2 = max(xv);
%         y1 = min(yv);
%         y2 = max(yv);
%         xbuff = [x1 x2 x2 x1 x1]/1000;
%         ybuff = [y1 y1 y2 y2 y1]/1000;
%         plot(xbuff,ybuff,'-k')
        
               
        end
        
        
        if isfield(BATH,'MAT')
        
        indx1 = BATH.MAT.z<0;
        indx2 = BATH.MAT.z>=0;
        
            
        scatter(BATH.MAT.x(indx1)/1000,BATH.MAT.y(indx1)/1000,5,BATH.MAT.z(indx1))          
        cb=colorbar;
        ylabel(cb,'z (m)')
        hold on
        plot(BATH.MAT.x(indx2)/1000,BATH.MAT.y(indx2)/1000,'.k')
                
        hold on
        indx = mark>0;
        plot(xe(indx)/1000,ye(indx)/1000,'k.')
%         x1 = min(xv);
%         x2 = max(xv);
%         y1 = min(yv);
%         y2 = max(yv);
%         xbuff = [x1 x2 x2 x1 x1]/1000;
%         ybuff = [y1 y1 y2 y2 y1]/1000;
%         plot(xbuff,ybuff,'-k')
        end
        
        ylabel('y UTM (km)')
        xlabel('x UTM (km)')
        
        title('NCOM model nesting')
        
        print -djpeg -r300 figure_NCOM_nesting
        
        close
        
        
        end
        
        
    end
