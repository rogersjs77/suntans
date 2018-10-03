%%  Lee Wave Data Maker
% This script:
% 1) Reads in u, w, T, eta, depth, and q data from a SUNTANS lee_wave simulation 
% with witch, sine, or arbitrary bathymetry. 
% 2) computes the form drag on the topography and energy flux from the 
% mean flowinto the wave.
% 3) Saves all data to a .mat file for plotting with separate scripts. 
% 4) Optionally produce plots of P and box fluxes at three points in run as
% well as flux as a function of time.

%% General Settings
EMPTY=999999; % code in SUNTANS data for uncomputed cells... do not change.
% datadir ='../../runs/3-9_h200_nonhs_witch/'; % directory of SUNTANS data
savename = 'lee_wave_data'; % name of saved .mat file
savedir = datadir; % directory to save data .mat file

makeplots = 0; % zero is off, 1 is on

%% Bathymetry settings
% toggle between a witch bottom or a sine bottom for form drag calculation
% sinebottom = 0;
% or for arbitrary bathymetry, use readdepth==1 and sinebottom==0 to read 
% in depth from depth.dat-voro. This uses cubic interpolation to get dhdx 
% at the voronoi points
readdepth = 0;

% Bathymetry dependent parameters to set for a given run
% sinebottom:
lambda = 1*6075;
K = 2*pi/lambda;
% witch bottom:
Lw =  1*5000; %half width of witch of agnesi

% if readdepth == 1
    depth=load([datadir,'/depth.dat-voro']);
    depth = depth(:,3);
    a = max(depth)-min(depth);
%     D = sum(depth)/length(depth); %average depth
    D = max(depth); %max depth
% else
%     a = 10; 
%     D = 0*8200 + 0*16000 + 0*2187 + 1*6911; %2172.2; %
% end

%% Physical Constants
rho0 = 1000;
grav = 9.81;
tday = 24*3600;
u0 = .88;
N = 0.007;

%% Loading grid and time data
c = load([datadir,'/cells.dat']);
points = load([datadir,'/points.dat']);
e = load([datadir,'/edges.dat']);
dz = load([datadir,'/vertspace.dat']);
Nkmax = getvalue([datadir,'/suntans.dat'],'Nkmax');
z = getz(dz); %depth
[Nc,cols] = size(c); %[number of cells, number of columns in cells.dat]

dt = getvalue([datadir,'/suntans.dat'],'dt');
nsteps = getvalue([datadir,'/suntans.dat'],'nsteps');
ntout = getvalue([datadir,'/suntans.dat'],'ntout');
gamma = getvalue([datadir,'/suntans.dat'],'gamma');
tend = nsteps*dt;

% Grid
xv = c(:,2);
yv = c(:,3);
[xvs,is] = sort(xv);
dx = xv(2)-xv(1);
dy = 2*yv(1);
Z = (z*ones(1,Nc))'; % z for every cell
ZZ = D+Z;
dZ = (dz*ones(1,Nc))';
dV = dZ.*dx.*dy; % cell volumes


% calculate domain length from grid data
Lx = max(xv)+dx/2;
Ly = max(yv)+dy/2;
xc = Lx/2; %witch location
X = (xv-xc)*ones(1,Nkmax);


%% Bathymetry calculations
if sinebottom==1
    h = 0.5*a.*cos(K.*xv);
    dhdx = -0.5*a*K.*sin(K.*xv);
elseif readdepth==1
    h = depth;
    dhdx = diff(h)/dx;
    xdiff = (xv(1:end-1)+dx/2);
    dhdx = interp1(xdiff, dhdx, xv, 'linear', 'extrap');
else
    % define the witch of agnesi
    h = a./(1+((xv-xc)/Lw).^2);
    % also will need dh/dx of witch
    dhdx = a/dx.*(1./(1+((xv-xc+dx/2)/Lw).^2) - 1./(1+((xv-xc-dx/2)/Lw).^2));
end

%% Flux box locations
%Horizontal locations of flux box
if sinebottom == 1 || readdepth==1
    i1 = 1;
    i2 = Nc;
else
    Lbox = 2*pi*Lw; % this puts the box sides around where h = 0.1*hmax
    x1 = xc-Lbox;
    x2 = xc+Lbox;
    i1 = find(abs(xv-x1)==min(abs(xv-x1)));
    i2 = find(abs(xv-x2)==min(abs(xv-x2)));
    i1 = min(i1);
    i2 = min(i2);
end

% Vertical location of flux box
h_flux = 1/4*2*pi*u0/N; %one quarter wavelength above floor, first maximum in pw_Gill
quarts0 = 1;
quarts = quarts0;
while h_flux < a % loop to ensure that the top of the box is above the maximum bathymetric height
    h_flux = h_flux + 1/16*2*pi*u0/N;
    quarts = quarts+1;
    disp('hill was larger than the initial box height')
    disp(['box height adjusted to be ',...
        num2str(quarts), '-quarter wavelengths tall (', num2str(h_flux),'m)'])
end
k_flux = min(find(abs(z-(-D+h_flux))==min(abs(z-(-D+h_flux)))));

%% Declare temperature and velocity data id (and compute nout)
% indexed as u(cell_id,level_id,timestep)
% this is all cell centered, so the e-flux calculations using it are less
% accurate than the full data due to interpolation of edge data to cell
% centers
Tfid=fopen([datadir,'/T.dat'],'rb');
nout = floor(length(fread(Tfid,'float64'))/(Nc*Nkmax)); 
t=[1:nout]*dt*ntout;
fclose(Tfid);

ufid=fopen([datadir,'/u.dat'],'rb');
Tfid=fopen([datadir,'/T.dat'],'rb');
eta_fid=fopen([datadir,'/fs.dat'],'rb');
qfid=fopen([datadir,'/q.dat'],'rb');


%% Linear Theory
% Bell's solution for a witch bottom (tidal flow)
Ly = dy;
h0 = a;

% Gill's solution (for both the sine wave and the witch bottom... the
% latter is actually Queney's)
if sinebottom == 1
   h0 = 0.5*a;
    %solution for a sine wave bottom (6.8 page 143 equation 6.8.7)
    pw_Gill = 1/2*K*rho0*h0^2*u0^2*sqrt(N^2-u0^2*K^2);
    E_Gill = pw_Gill*Lx*Ly;
    %field solutions
    M = sqrt(N^2/u0^2 - K^2);
    p_gill = -u0^2*rho0*h0*(N^2/u0^2 - K^2)^(1/2).*sin(K.*X + (N^2/u0^2 - K^2)^(1/2).*ZZ);
    % from Khatiwala
    p_rigid = -rho0*u0^2*h0*M/sin(M*D).*cos(M.*(D-ZZ)).*cos(K.*X);
    w_gill = -u0*K*h0.*sin(K.*X + (N^2/u0^2 - K^2)^(1/2).*ZZ);
    pw_gill = p_gill .* w_gill;
    
else
    %solution for a witch bottom
    E_Gill = pi/4 * rho0 * u0^2 * Ly * N * h0^2;
    pw_Gill = rho0*u0^2*h0^2*N/Lw;
    %field solutions
    p_queney = -rho0*u0*N*a*Lw .* (Lw^2 + X.^2).^(-1) .* (Lw.*sin(N/u0.*ZZ) +...
    X.*cos(N/u0.*ZZ));
    w_queney = u0*a*Lw.* (Lw^2 + X.^2).^(-2) .*...
    [sin(N/u0.*ZZ).*(X.^2-Lw^2) - 2*Lw.*cos(N/u0.*ZZ).*X];
    pw_queney = p_queney .* w_queney;
end

%% Initialize arrays to fill
% energy arrays to fill
Ek = zeros(nout,1);
Ep = Ek;
Etot = Ek;
Ekprime = Ek;
Epprime = Ek;

% flux box arrays to fill
F_t = zeros(nout,1);
F_r = zeros(nout,1);
F_l = zeros(nout,1);

% form drag related arrays to fill
force_per_y = zeros(nout,1);
power_per_y = zeros(nout,1);
cell_center_power = zeros(nout,1);

% array of pw(z) at x=xc, where pw is the correlation of perturbation
% pressure and vertical velocity (so the energy flux per area)
pw_on_hill = zeros(Nkmax,nout);

% arrays of p, w, u, and T for field animations... warning, these can be
% large
panimate = zeros(ntout,Nc,Nkmax);
wanimate = panimate;
tanimate = panimate;

%% The Heavy Lifting:
% Read in data for each nout
% compute form drag and energy fluxes
% fill the above arrays
for n=1:nout
    udata=fread(ufid,Nc*Nkmax*3,'float64');
    ud = reshape(udata,Nc,Nkmax,3);
    u = squeeze(ud(:,:,1));
    w = squeeze(ud(:,:,3));
    Tdata=fread(Tfid,Nc*Nkmax,'float64');
    T = reshape(Tdata,Nc,Nkmax);
    eta = fread(eta_fid,Nc,'float64');
    
    % for the first nout, define bottom cells, fix dZ in bottom cells
    % and the background temperature profile, T_bar
    if(n==1)
        empties = zeros(size(u));
        empties(u==EMPTY)=1;
        u(u==EMPTY)=nan;
        w(w==EMPTY)=nan;
        T(T==EMPTY)=nan;
        bottom_cells = zeros(Nc,1);
            for i = 1:Nc
                if isnan(T(i,Nkmax))
                    bottom_cells(i) = find(isnan(T(i,:)),1)-1;
                else
                    bottom_cells(i) = Nkmax;
                end
                bot = bottom_cells(i);
                dZ(i,bot) = depth(i)+(z(bot)+dz(bot)/2);
            end  
        % fix any empty cells now to 0
        u(isnan(u))=0;
        w(isnan(w))=0;
        T(isnan(T))=0;
        
        % save this w for animation
        wanimate(n,:,:) = w;

        %compute initial temperature profile
        T0 = 20 + (N^2/(gamma*grav)).*z'; %computed using N, as was done in SUNTANS
        Tdiff = T(1,:)-T0; %compare initial data T to this theoretical setting (should be spot on)
        %     figure();
        %     plot(Tdiff,z)
        teval = find(depth == max(depth),1); % make sure the T_bar taken from the data is from the deepest place!
        T_bar=ones(Nc,1)*T(teval,:);
        %T_bar = ones(Nc,1)*T0;
    end
    
    
    % if n ~=1, still have to fix any empty cells to 0's
    u(u==EMPTY)=0;
    w(w==EMPTY)=0;
    T(T==EMPTY)=0;
    wanimate(n,:,:) = w;
    
    % density
    rho = rho0.*(-gamma.*T);
    rho_prime = -gamma*rho0*(T-T_bar);
    
    % hydrostatic pressure (with/without use of free surface
    % deflections and/or cut cell bottm dz's)
    p_prime = zeros(Nc,Nkmax);
    p_prime(:,1) = 0.5*rho_prime(:,1)*grav.*(dz(1)+2*eta.*1e8); % including the free surface in the pressure
%     p_prime(:,1) = 0.5*rho_prime(:,1)*grav.*(dz(1));
%     for k=2:Nkmax
%         p_prime(:,k) = p_prime(:,k-1)+0.5*rho_prime(:,k-1)*grav*dz(k-1)+...
%                                       0.5*rho_prime(:,k)*grav*dz(k);
%     end
%    here's a p_prime calculation using cut cell dz values at the bottom
    for i=1:Nc
        for k = 2:bottom_cells(i)
        p_prime(i,k) = p_prime(i,k-1)+0.5*rho_prime(i,k-1)*grav*dZ(i,k-1)+...
                                      0.5*rho_prime(i,k)*grav*dZ(i,k);
        end
    end
    
    % nonhydrostatic pressure
    q = reshape(fread(qfid,Nc*Nkmax,'float64'),Nc,Nkmax);
    % add nonhydrostatic pressure to p_prime
    p_prime = p_prime+q.*rho0;
    % clear out any pressure in non-computational cells
    p_prime(empties==1) = 0;
    
    % fill animation arrays
    panimate(n,:,:) = p_prime;
    tanimate(n,:,:) = T;

    % Barotropic currents
    U = sum(u.*dZ,2)./sum(dZ,2)*ones(1,Nkmax);
     
    % compute p_prime at the floor (rather than cell center in bottom cell) 
    % and calculate form drag related quantities 
    for i = 1:Nc
        kb = bottom_cells(i);
%         dZ(i,kb) = dZ(i,kb) - h(i) + (D+z(kb)-0.5*dz(kb));
        
        % p_prime_bottom = p_prime(i,kb-1) + rho_prime(i,kb)*dZ(i,kb); 
        p_prime_bottom = p_prime(i,kb) + rho_prime(i,kb)*0.5*dZ(i,kb); 
        
        center_force_per_area = p_prime(i,kb)*dhdx(i);
        force_per_area = p_prime_bottom*dhdx(i);
        power_per_area = U(i,kb)*force_per_area;
        force_per_y(n) = force_per_y(n) + force_per_area * dx;
        power_per_y(n) = power_per_y(n) + power_per_area * dx;
        cell_center_power(n) = cell_center_power(n) + U(i,kb)*center_force_per_area*dx*dy;
    end

    pw = p_prime.*w.*dx*dy;
    pu = p_prime.*(u-U).*dZ.*dy;
    Fv = sum(pw,1);
    Fh = sum(pu,2);

    % Box flux
    F_l(n) = sum(pu(i1,k_flux:Nkmax));
    F_r(n) = sum(pu(i2,k_flux:Nkmax));
    F_t(n) = sum(pw(i1:i2,k_flux));    
    
    %  add to pw_on_hill
    pw_on_hill(:,n) = squeeze(p_prime(floor(Nc/2),:).*w(floor(Nc/2),:));

    % energy calculations
    Ek(n) = sum(sum((u.*u + w.*w).*rho.*dV));
    Ep(n) = sum(sum(grav.*rho.*(D+Z).*dV));
    Etot(n) = Ek(n)+Ep(n);
    Ekprime(n) = sum(sum(((u-U).*(u-U) + w.*w).*rho.*dV));
    Epprime(n) = sum(sum(grav.*rho_prime.*(D+Z).*dV));
    
    
    % plotting p' and fluxes at three different times in run
    if (makeplots==1 && (n== floor(nout/3) || n== floor(nout/2) || n==nout))   
        figure()
        subplot(3,2,[1,2])
        contourf(X,Z,p_prime,'LineStyle','none')
        colorbar
        title(['perturbation pressure (kg/m/s^2) at t=',num2str(n*ntout*dt/tday,3),' days'])
        xlabel('easting (m)')
        ylabel('depth (m)')
        hold on
        plot((xv(i1)-xc).*ones(size(k_flux:Nkmax)),z(k_flux:Nkmax),'r',...
            xv(i1:i2)-xc,z(k_flux).*ones(size(i1:i2)),'r',...
            (xv(i2)-xc).*ones(size(k_flux:Nkmax)),z(k_flux:Nkmax),'r')
        hold off
        
        subplot(3,2,[3,4])
        plot(xv(i1:i2)-xc,pw(i1:i2,k_flux)./(pw_Gill*dx*dy))
        title(['p*w/pw_{Gill} along top of box, z=',num2str(D+z(k_flux)),'m'])
        xlabel('easting (m)')
        ylabel('p*w/pw_{Gill}')
        
        subplot(3,2,5)
        plot(pu(i1,k_flux:Nkmax)./(dZ(k_flux:Nkmax).*dy*pw_Gill),D+z(k_flux:Nkmax))
        title(['p*u along left of box, x=',num2str(xv(i1)-xc),'m'])
        xlabel('p*u/pw_{Gill}')
        ylabel('depth (m)')
        
        subplot(3,2,6)
        plot(pu(i2,k_flux:Nkmax)./(dZ(k_flux:Nkmax).*dy*pw_Gill),D+z(k_flux:Nkmax))
        title(['p*u along right of box, x=',num2str(xv(i2)-xc),'m'])
        xlabel('p*u/pw_{Gill}')
        ylabel('depth (m)')
        
    end
    
end

force_from_hill = force_per_y.*dy;
power_from_hill = power_per_y.*dy;

%% Saving
save([savedir,savename])

%% Plotting energy flux vs time
if makeplots==1
    figure()  
    subplot(2,1,1)
    plot(t/tday,F_r/E_Gill,'k:',...
         t/tday,F_l/E_Gill,'k-.',...
         t/tday,F_t/E_Gill,'k--',...
         t/tday,(F_r-F_l+F_t)/E_Gill,'k-');
    title('Fluxes out of a box around the hill');
    legend(sprintf('Right x=+%.1f km',(xv(i1)-xc)/1000),...
           sprintf('Left x=%.1f km',(xv(i2)-xc)/1000),...
           sprintf('Top z=%.1f m',h_flux),...
           'Net (Top + Right - Left)',...
           'location','east');
    xlabel('Time (d)');
    ylabel('Flux/F_{Gill} @ z_{flux}');

    subplot(2,1,2)
    plot(t/tday,power_from_hill/E_Gill,'k.-',...
        t/tday,cell_center_power/E_Gill,'k:',...
        t/tday,ones(size(t)),'k--');
    title('Flux via form drag');
    xlabel('Time (d)');
    ylabel('Flux/F_{Gill} @ z_{flux}');
    legend('p at cell bottom','p at cell center','location','east')
end

%% Cleaning up
fclose all;

