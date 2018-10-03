function [] = quadgrid_periodic(datadir,L,W,Nx,Ny,BC,STRETCHING,CHEBYCHEV,Lr,K,rmax,x0,y0,theta)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% File name: quadrid.m
% Description: Create Cartesian grid files for use with the quad
% version of suntans.  Horizontal grid stretching can be used by
% setting STRETCH=true;
%
% If stretching in the x-direction is employed, then the grid consists
% of a refined region in the middle of the domain with a width Lr and
% number of grid points Nxr in this region.  Nxr can be set either by
% deciding on a desired resolution in the refined region or it can
% just be set arbitrarily. i.e. There is no need to set dxr (see
% below) to compute a stretched grid.
%
% The grid is stretched to the left and right of the refined 
% region by an amount r which is determined to ensure that the
% total length of the domain matches L.  Schematically, the grid
% centers would look like
%
%           | .    .   .  . ........... .  .   .    . |
%
%           |<-----Ls ----->|<---Lr-->|<------Ls ---->|
%
% Here, the length of the refined region is Lr and the length of each
% of the stretched regions is Ls=(L-Lr)/2.  Therefore, the number of
% grid cells in the refined region is Nxr and the number in each
% stretched region is Nxs=(Nx-Nxr)/2.  Note that if (Nx-Nxr) is not 
% divisible by 2 then Nxr is increased by 1.
%
% The function stretchgrid(xpg,ypg,Lr,Nxr,rmax) takes as its
% input xpg and ypg which are the vertices of the grid 
% arranged in 2D arrays (i.e. not 1d arrays). rmax should be 1.1
% but it can of course be large depending on how mutch stretching
% is desired. Note that in order for r to be obtained, a solution
% to the following algebraic equation is needed:
%
% r^Nxs - 1 - Ls*(Nxr/Lr)*(r-1) = 0,
%
% In cases of small Nxs (too few cells in stretched region) or
% small Lr/Nxr (too much refinement), this may require an
% exceedingly large stretching factor and in some cases the solver
% fsolve() will not find a solution.
%
% revised by Yun Zhang 2/18/2013 @Stanford
% 1) make dx=dy=1 and dxr=1/K to avoid numerical error fault
% 2) make sure grad(i,1)!=-1
%
% revised by Justin Rogers 5/2018 @Stanfod
% allow for periodic boundary conditions
% allow for x0 = [x0,y0], lower left corner coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Directory in which points.dat, cells,dat, edges.dat files will be
% placed
% datadir='./'

% Length and width of domain
% L = 25000;
% W = 1000;

% Number of cells
% Nx = 100;
% Ny = 1;
delete(gcp('nocreate'))
poolobj = parpool('local'); 

ampfacx=Nx/L;
ampfacy=Ny/W;
L=Nx;
W=Ny;
dx=1;
dy=1;

Mixed = BC>5; % find mixed type2/3 BC
BC(Mixed)=2; % set BC back to type 2

if nargin<12
   x0=0; % lower left corner x0=E, 
   y0=0; % y0=N;
   theta=0;
end
if nargin<14
    theta=0;
end
% Whether or not to employ stretching
% STRETCHING=false;
% 
% % If stretching is employed, use Chebychev grids
% CHEBYCHEV=false;

% Length of resolved region
Lr_x = Lr(1);
if length(Lr)>1
    refine2D=1;
    Lr_y = Lr(2);
else
    refine2D=0;
    Lr_y=0;
end

Lramp_x =Lr_x*ampfacx;
Lramp_y =Lr_y*ampfacy;
% Resolution of resolved region (which is constant). In this case
% the grid spacing is half the original spacing.
Kx=K(1); % the resolution for dx/dxr
if refine2D
    Ky=K(2);  % the resolution for dy/dyr
else
    Ky=1;
end
dxr = 1/Kx;
dyr = 1/Ky;

% Number of grid points in resolved region
Nxr = ceil(Lramp_x/dxr);
Nyr = ceil(Lramp_y/dyr);

% Maximum acceptable stretching factor. The code will exit if the
% stretching factor needed to give the refined region exceeds rmax.
% rmax = 1.1;

% Boundary condition types:
% 1 solid free-slip wall
% 2 velocity specified
% 3 free-surface specified
% 5 periodic (converted back to 0 afterwards)
EastBoundary = BC(1);
NorthBoundary = BC(2);
WestBoundary = BC(3);
SouthBoundary = BC(4);
% WestBoundary = 5;
% EastBoundary = 5;
% NorthBoundary = 1;
% SouthBoundary = 1;

N = (Nx+1)*(Ny+1);
[xpg,ypg] = ndgrid([0:dx:L],[0:dy:W]);
xp = xpg(:);
yp = ypg(:);
Np = length(xp);

[xv,yv] = ndgrid([dx/2:dx:L-dx/2],[dy/2:dy:W-dy/2]);
xv = xv(:);
yv = yv(:);
Nc = length(xv);
 
cells = zeros(Nc,4);

parfor n=1:Nc
  xcell = xv(n) + [-dx/2,-dx/2,dx/2,dx/2];
  ycell = yv(n) + [-dy/2,dy/2,dy/2,-dy/2];
  for m=1:4
    cells(n,m) = find(xp==xcell(m) & yp==ycell(m));
  end
end  

[xeu,yeu] = ndgrid([0:dx:L],[dy/2:dy:W-dy/2]);
[xev,yev] = ndgrid([dx/2:dx:L-dx/2],[0:dy:W]);
xeu = xeu(:);
yeu = yeu(:);
xev = xev(:);
yev = yev(:);

Ne = Nx*(Ny+1) + Ny*(Nx+1);
mark = zeros(Ne,1);
edges = zeros(Ne,2);
grad = -ones(Ne,2);

k=1;
for n=1:length(xeu)
  edges(k,1) = find(xp==xeu(n) & yp==yeu(n)-0.5*dy);
  edges(k,2) = find(xp==xeu(n) & yp==yeu(n)+0.5*dy);
  if(xeu(n)==0)
    mark(k)=WestBoundary;
  elseif(xeu(n)==L)
    mark(k)=EastBoundary;
  else
    mark(k)=0;
  end

  xv1 = xeu(n)-0.5*dx;
  xv2 = xeu(n)+0.5*dx;
  yv1 = yeu(n);
  yv2 = yeu(n);

  if(WestBoundary==5 & EastBoundary==5)
    if(xv1<0)
      xv1=L-0.5*dx;
    end	
    if(xv2>L)
      xv2=0.5*dx;
    end
  end
  
  ind1 = find(xv1==xv & yv1==yv);
  ind2 = find(xv2==xv & yv2==yv);
  if(~isempty(ind1))
    grad(k,1) = ind1;
  end
  if(~isempty(ind2))
    grad(k,2) = ind2;
  end

  if(isempty(ind1) & isempty(ind2))
    error('Couldn''t find at least one neighbor for edge %d!\n',k);
  end

  k=k+1;
end
for n=1:length(xev)
  edges(k,1) = find(xp==xev(n)-0.5*dx & yp==yev(n));
  edges(k,2) = find(xp==xev(n)+0.5*dx & yp==yev(n));
  if(yev(n)==0)
    mark(k)=SouthBoundary;
  elseif(yev(n)==W)
    mark(k)=NorthBoundary;
  else
    mark(k)=0;
  end

  xv1 = xev(n);
  xv2 = xev(n);
  yv1 = yev(n)+0.5*dy;
  yv2 = yev(n)-0.5*dy;

  if(NorthBoundary==5 & SouthBoundary==5)
    if(yv1>W)
      yv1=0.5*dy;
    end	
    if(yv2<0)
      yv2=W-0.5*dy;
    end
  end

  ind1 = find(xv1==xv & yv1==yv);
  ind2 = find(xv2==xv & yv2==yv);
  if(~isempty(ind1))
    grad(k,1) = ind1;
  end
  if(~isempty(ind2))
    grad(k,2) = ind2;
  end

  if(isempty(ind1) & isempty(ind2))
    error('Couldn''t find at least one neighbor for edge %d!\n',k);
  end

  k=k+1;
end

for j=1:Ne
  if(grad(j,1)==-1 | grad(j,2)==-1)
    if(grad(j,1)==-1)
      grad(j,1)=grad(j,2);
      grad(j,2)=-1;
    end
  end
end  

neigh = -ones(Nc,4);
for n=1:Nc
  % clockwise get neigh to make we can get share node for each boundary
  ycell = yv(n);
  xcell = xv(n) - dx;
  if(WestBoundary==5 & xcell<0)
    xcell = L-dx/2;
  end

  ind = find(xcell==xv & ycell==yv);
  if(~isempty(ind))
    neigh(n,1)=ind;
  end

  xcell = xv(n);
  ycell = yv(n) + dy;
  if(NorthBoundary==5 & ycell>W)
    ycell = dy/2;
  end

  ind = find(xcell==xv & ycell==yv);
  if(~isempty(ind))
    neigh(n,2)=ind;
  end

  ycell = yv(n);
  xcell = xv(n) + dx;
  if(EastBoundary==5 & xcell>L)
    xcell = dx/2;
  end

  ind = find(xcell==xv & ycell==yv);
  if(~isempty(ind))
    neigh(n,3)=ind;
  end

  xcell = xv(n);
  ycell = yv(n) - dy;
  if(SouthBoundary==5 & ycell<0)
    ycell = W-dy/2;
  end

  ind = find(xcell==xv & ycell==yv);
  if(~isempty(ind))
    neigh(n,4)=ind;
  end
end

%% now set mixed type2/3 bc

% set cell =type3 on left looking into domain
%E,N,W,S
if max(Mixed)>0
    xe = [xeu; xev];
    ye = [yeu; yev];    
    rmin = 2; % radius 2 cells
    if Mixed(1) % E = lower right
        r = sqrt((xe-max(xe)).^2+(ye-min(ye)).^2);
        indx = r<rmin & mark==2; % within radius and type2
        mark(indx)=3;
    end
    if Mixed(2) %N=upper right
        r = sqrt((xe-max(xe)).^2+(ye-max(ye)).^2);
        indx = r<rmin & mark==2; % within radius and type2
        mark(indx)=3;
    end
    if Mixed(3) % W upper left
        r = sqrt((xe-min(xe)).^2+(ye-max(ye)).^2);
        indx = r<rmin & mark==2; % within radius and type2
        mark(indx)=3;
    end
    if Mixed(4) % S lower left
        r = sqrt((xe-min(xe)).^2+(ye-min(ye)).^2);
        indx = r<rmin & mark==2; % within radius and type2
        mark(indx)=3;
    end
    
end

% 
% figure(1);
% clf;
% hold on;
% axis([-dx/2 L+dx/2 -dy/2 W+dy/2]);
% axis off;
% 
% for n=1:Ne
%   if(mark(n)==0)
%     color='k.-';
%   else
%     color='r.-';
%   end
%   
%   plot(xp(edges(n,:))',yp(edges(n,:))',color);
%   if(grad(n,1)~=-1 & grad(n,2)~=-1)
%     plot(xv(grad(n,:)),yv(grad(n,:)),'m-');
%   end
% end
%   
% for n=1:Nc
%   for m=1:4
%     if(neigh(n,m)~=-1)
%       plot([xv(n),xv(neigh(n,m))],[yv(n),yv(neigh(n,m))],'b.-');
%     end
%   end
% end

% Now we can stretch the grid without destroying the connectivity
if(STRETCHING)
  if(~CHEBYCHEV)
    [xv,yv,xp,yp]=stretchgrid(xpg,ypg,Lramp_x,Nxr,rmax);
    if refine2D % refinement in y direction too
        [yvt,xvt,ypt,xpt]=stretchgrid(reshape(yp,Nx+1,Ny+1)',reshape(xp,Nx+1,Ny+1)',Lramp_y,Nyr,rmax);
        % put back in same order for connectivity
        xv = reshape(xvt,Ny,Nx)';
        xv = xv(:); 
        yv = reshape(yvt,Ny,Nx)';
        yv = yv(:);
        xp = reshape(xpt,Ny+1,Nx+1)';
        xp = xp(:);        
        yp = reshape(ypt,Ny+1,Nx+1)';
        yp = yp(:);
%         xv=xvt;
%         yv=yvt;
%         xp=xpt;
%         yp=ypt;
    end
  else
    j_cheb = [1/2:Nx+1/2];
    x_cheb = 0.5*(1-cos((j_cheb-1/2)*pi/Nx));

    [xp,yp] = ndgrid(x_cheb,[0:1/Ny:1]);
    fprintf('dx min = %f\n',min(diff(x_cheb)));
    j_cheb = [1:Nx];
    x_cheb = 0.5*(1-cos((j_cheb-1/2)*pi/Nx));

    [xv,yv] = ndgrid(x_cheb,[0.5:1/Ny:1-0.5]);

    xp = ampfacx*xp(:);
    yp = ampfacy*yp(:);
    xv = ampfacx*xv(:);
    yv = ampfacy*yv(:);
  end
end

% Remove repeated edges on periodic boundaries
xe=mean(xp(edges(:,[1,2])),2)/ampfacx;
ye=mean(yp(edges(:,[1,2])),2)/ampfacy;
% rotate to absolute ref frame
xe2 = xe*cosd(theta)-ye*sind(theta);
ye2 = xe*sind(theta)+ye*cosd(theta);
% locate lower left corner at x0,y0
xe = xe2 + x0;
ye = ye2 + y0;
clear xe2 ye2

% xe=mean(xp(edges(:,[1,2])),2)+x0;
% ye=mean(yp(edges(:,[1,2])),2)+y0;
inds = find((mark==5 & xe==L) | (mark==5 & ye==W));
grad(inds,:) = [];
edges(inds,:) = [];
mark(inds,:) = [];
Ne = length(edges(:,1));

% Set Periodic edges to type 0 as they are computational edges
mark(find(mark==5))=0;

% veronoi output
xv_n = xv;
yv_n = yv;
xv = xv./ampfacx;
yv = yv./ampfacy;
% rotate to absolute ref frame
xv2 = xv*cosd(theta)-yv*sind(theta);
yv2 = xv*sind(theta)+yv*cosd(theta);
% locate lower left corner at x0,y0
xv = xv2 + x0;
yv = yv2 + y0;
clear xv2 yv2

celloutput=zeros(Nc,11);
celloutput(:,1)=4;
celloutput(:,2)=xv;
celloutput(:,3)=yv;
celloutput(:,4:7)=cells;
celloutput(:,8:11)=neigh;

edgeoutput=zeros(Ne,5);
%find open BC to make sure grad[2*j]~=-1
loc=find(mark~=0 & mark~=1);
edgeoutput(:,1:2)=edges;
edgeoutput(:,3)=mark;
edgeoutput(:,4:5)=grad;
for m=1:length(loc)
    if grad(m,1)==-1
        edgeoutput(m,5)=-1;
        edgeoutput(m,4)=grad(m,2);
    end
end

% points
xp = xp/ampfacx;
yp = yp/ampfacy;
% rotate to absolute ref frame
xp2 = xp*cosd(theta)-yp*sind(theta);
yp2 = xp*sind(theta)+yp*cosd(theta);
% locate lower left corner at x0,y0
xp = xp2 + x0;
yp = yp2 + y0;
clear xp2 yp2 

pointoutput=zeros(Np,3);
pointoutput(:,1)=xp;
pointoutput(:,2)=yp;

%check results (should be edge length)
dis=3*ones(Ne,1);
for i=1:Ne
    nc1=grad(i,1);
    nc2=grad(i,2);
    if(nc1~=-1 & nc2~=-1)
        dis(i)=((xv_n(nc1)-xv_n(nc2))^2+(yv_n(nc1)-yv_n(nc2))^2)^0.5;
    end
end



celloutput(:,4:7)=celloutput(:,4:7)-1;
for i=8:11
loc=find(celloutput(:,i)~=-1);
celloutput(loc,i)=celloutput(loc,i)-1;
end

edgeoutput(:,1:2)=edgeoutput(:,1:2)-1;

for i=4:5
loc=find(edgeoutput(:,i)~=-1);
edgeoutput(loc,i)=edgeoutput(loc,i)-1;
end


cells_file = [datadir,'/cells.dat'];
cellsf = fopen(cells_file,'w');
fprintf(cellsf, '%1.0f %12.10e %12.10e %1.0f %1.0f %1.0f %1.0f %1.0f %1.0f %1.0f %1.0f\n', celloutput');
status = fclose(cellsf);   

points_file = [datadir,'/points.dat'];
pointsf = fopen(points_file,'w');
fprintf(pointsf, '%12.10e %12.10e %1.0f\n', pointoutput');
status = fclose(pointsf);   

edges_file = [datadir,'/edges.dat'];
edgesf = fopen(edges_file,'w');
fprintf(edgesf, '%1.0f %1.0f %1.0f %1.0f %1.0f\n', edgeoutput');
status = fclose(edgesf);

save SUNTANS_grid_quadgrid -v7.3

delete(poolobj)

