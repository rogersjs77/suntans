function [] = SUNTANS_process(datadir,istep,Tavg,quad)
%% SUNTANS Process
% This script reads standard SUNTANS output and saves to a matlab .mat file
% for further processing
%
% Justin Rogers
% Stanford University
% July 2017
%
%
%% General Settings
disp('processing data')
% clear
% datadir='../data';
% state = 'sal'; %'temp','sal','full'
EMPTY=999999; % code in SUNTANS data for uncomputed cells... do not change.
% datadir ='../../runs/3-9_h200_nonhs_witch/'; % directory of SUNTANS data
% savename = 'SUNTANS_results'; % name of saved .mat file
savedir = './'; % directory to save data .mat file

% makeplots = 0; % zero is off, 1 is on
rho0 = 1000;
grav = 9.81;

%% Bathymetry settings

depth=load([datadir,'/depth.dat-voro']);
xv = depth(:,1);
yv = depth(:,2);
depth = depth(:,3);
D = max(depth); %max depth


%% Loading grid and time data
c = load([datadir,'/cells.dat']);
points = load([datadir,'/points.dat']);
edges = load([datadir,'/edges.dat']);
dz = load([datadir,'/vertspace.dat']);
Nkmax = getvalue([datadir,'/suntans.dat'],'Nkmax');
z = getz(dz); %depth
[Nc,~] = size(c); %[number of cells, number of columns in cells.dat

dt = getvalue([datadir,'/suntans.dat'],'dt');
nsteps = getvalue([datadir,'/suntans.dat'],'nsteps');
ntout = getvalue([datadir,'/suntans.dat'],'ntout');
gamma = getvalue([datadir,'/suntans.dat'],'gamma');
beta = getvalue([datadir,'/suntans.dat'],'beta');
sponge_distance = getvalue([datadir,'/suntans.dat'],'sponge_distance');
CdB = getvalue([datadir,'/suntans.dat'],'CdB');
z0B = getvalue([datadir,'/suntans.dat'],'z0B');
nu_H = getvalue([datadir,'/suntans.dat'],'nu_H');
nu = getvalue([datadir,'/suntans.dat'],'nu');
thetaramp = getvalue([datadir,'/suntans.dat'],'thetaramptime');
basetime = getvalue([datadir,'/suntans.dat'],'basetime');
starttime = getvalue([datadir,'/suntans.dat'],'starttime');
wave_nesting = getvalue([datadir,'/suntans.dat'],'wave_nesting');

tend = nsteps*dt;

[xvs,is] = sort(xv);
dx = xv(2)-xv(1);
dy = 2*yv(1);
Z = (z*ones(1,Nc))'; % z for every cell
ZZ = D+Z;
dZ = (dz*ones(1,Nc))';
dV = dZ.*dx.*dy; % cell volumes

%% find cut cell bottom
zb = z-dz/2;
zt = z+dz/2;
dZ_cut = nan+dZ;

for i=1:Nc
        % cells below depth
%         dZ_cut(i,:) = 0*z;
        % cells fully open
        indx = zb>=-depth(i);
        dZ_cut(i,indx)= zt(indx)-zb(indx);
        % cells partially closed
        indx = zb<-depth(i) & zt>-depth(i);
        dZ_cut(i,indx)= zt(indx)+depth(i);
end


%%


% calculate domain length from grid data
Lx = (max(xv)-min(xv))+dx/2;
Ly = (max(yv)-min(yv))+dy/2;
xc = Lx/2; %witch location
X = (xv-xc)*ones(1,Nkmax);

% get sponge and wavemaker locations
x_sponge = [sponge_distance, Lx-sponge_distance];

% this is the old style
f2 = importdata('../sources.c');
iforce = strcmp('  REAL dW; ',f2);
if max(iforce)==1
dW = [];    j=1; 
    while isempty(dW) && j<length(f2)
        [dW] = sscanf(f2{j},'  dW = %e;');
        j=j+1;
    end
else % now default
    dW = getvalue([datadir,'/suntans.dat'],'dW');
    if isempty(dW)
        dW=0;
    end
%     dW = sponge_distance;
end
clear fw iforce j
x_wavemaker = [x_sponge(1)+dW, x_sponge(2)-dW];

%% load in results
% figure out if netcdf (nc) or binary file format

nc_files = dir([datadir '/sun_out.nc*']);
nc_AVG = dir([datadir '/sun_AVG.nc*']);
binary_files = dir([datadir '/u.dat*']);

if ~isempty(nc_files)
    DATA = netcdf_load_struct([datadir '/' nc_files.name]);
    t_ref = DATA.time;
    t = t_ref-min(t_ref);
    Nout = length(t_ref);
    t = t(1:istep:Nout);
    mtime = t/24/3600+datenum(sprintf('%.6f',starttime),'yyyymmdd.HHMMSS');
    Nout_keep = length(t);
    
    step=1:istep:Nout;
    u = DATA.uc(:,:,step);
    v = DATA.vc(:,:,step);
    w = DATA.w(:,:,step);
    T = DATA.temp(:,:,step);
    S = DATA.salt(:,:,step);
    q = DATA.q(:,:,step);
    nuT = DATA.nu_v(:,:,step);
    Uface = DATA.U(:,:,step);   
     
    eta(1:Nc,1,1:Nout_keep) = DATA.eta(:,step);
    
    xe = DATA.xe;
    ye = DATA.ye;
    n1 = DATA.n1;
    n2 = DATA.n2;
    
    clear DATA
    
    % load in Averages file    
    if ~isempty(nc_AVG)
        DATA = netcdf_load_struct([datadir '/' nc_AVG.name]);
        
        
        AVG.t_ref = DATA.time;
        AVG.t = AVG.t_ref - min(AVG.t_ref);
        AVG.mtime = AVG.t/24/3600+datenum(sprintf('%.6f',starttime),'yyyymmdd.HHMMSS');
        % take time rounding same as istep
        if length(AVG.t)>3
            AVG.istep = ceil((istep*dt*ntout)./(AVG.t(2)-AVG.t(1)));
        else
            AVG.istep=1;
        end
        AVG.Ntout = length(AVG.t_ref);
        step=1:AVG.istep:AVG.Ntout;
        if length(step)<2
            step = [1 AVG.Ntout];
        end
        
        AVG.mtime = AVG.mtime(step);
        AVG.t = AVG.t(step);
        AVG.t_ref = AVG.t_ref(step);
        AVG.Ntout = length(AVG.t_ref);
        
        AVG.u = DATA.uc(:,:,step);
        AVG.v = DATA.vc(:,:,step);
        AVG.w = DATA.w(:,:,step);
        AVG.nuT = DATA.nu_v(:,:,step);
        AVG.S = DATA.salt(:,:,step);
        AVG.T = DATA.temp(:,:,step);
%         AVG.q = DATA.q(:,:,step);
        AVG.eta(1:Nc,1,1:AVG.Ntout) = DATA.eta_avg(:,step);
%         AVG.alphaw = DATA.alphaw;
        
        % edge values
%         for i=1:length(z)
%             for j=1:length(DATA.time)
%                 AVG.u_avg(:,i,j) = DATA.u_avg(:,i,j).*n1;
%                 AVG.v_avg(:,i,j) = DATA.u_avg(:,i,j).*n2;
%                 
%                 AVG.u_avg(n1==0,i,j)=nan;
%                 AVG.v_avg(n2==0,i,j)=nan;
%                 
% %                 AVG.uw_var(:,i,j) = DATA.uw_var(:,i,j);
% %                 AVG.uw_var_avg(:,i,j) = abs(DATA.uw_var_avg(:,i,j).*n1);
% %                 AVG.vw_var_avg(:,i,j) = abs(DATA.uw_var_avg(:,i,j).*n2);
% %                 
% %                 AVG.uw_var_avg(n1==0,i,j)=nan;
% %                 AVG.vw_var_avg(n2==0,i,j)=nan;
%             end
%         end
        
        clear DATA
    end
    
elseif ~isempty(binary_files)

    % Declare temperature and velocity data id (and compute nout)
    % indexed as u(cell_id,level_id,timestep)
    % this is all cell centered, so the e-flux calculations using it are less
    % accurate than the full data due to interpolation of edge data to cell
    % centers
    Tfid=fopen([datadir,'/T.dat'],'rb');
    Nout = floor(length(fread(Tfid,'float64'))/(Nc*Nkmax)); 
    t=[1:Nout]*dt*ntout;
    t = t(1:istep:Nout);
    Nout_keep = length(t);
    fclose(Tfid);

    sfid=fopen([datadir,'/s.dat'],'rb');
    ufid=fopen([datadir,'/u.dat'],'rb');
    Tfid=fopen([datadir,'/T.dat'],'rb');
    etafid=fopen([datadir,'/fs.dat'],'rb');
    qfid=fopen([datadir,'/q.dat'],'rb');
    nufid=fopen([datadir,'/nut.dat'],'rb');

    % Load in results
    % Read in data for each nout
    precision='float64';
    % istep = 1;
    kk=1;
    for step=1:istep:Nout
        % velocity
        arraysize=Nc*Nkmax*3;
        data = getcdata_rev(ufid, arraysize, step, precision);
        data = reshape(data,Nc,Nkmax,3);
        u(:,:,kk) = squeeze(data(:,:,1));
        v(:,:,kk) = squeeze(data(:,:,2));
        w(:,:,kk) = squeeze(data(:,:,3));

        % temperature
        arraysize=Nc*Nkmax;
        data = getcdata_rev(Tfid, arraysize, step, precision);
        data = reshape(data,Nc,Nkmax);
        T(:,:,kk) = data;

        % eta
        arraysize=Nc;
        data = getcdata_rev(etafid, arraysize, step, precision);
        data = reshape(data,Nc,1);
        eta(:,:,kk) = data;

        % salinity
        arraysize=Nc*Nkmax;
        data = getcdata_rev(sfid, arraysize, step, precision);
        data = reshape(data,Nc,Nkmax);
        S(:,:,kk) = data;  

        % pressure
        arraysize=Nc*Nkmax;
        data = getcdata_rev(qfid, arraysize, step, precision);
        data = reshape(data,Nc,Nkmax);
        q(:,:,kk) = data;

        % nuT
        arraysize=Nc*Nkmax;
        data = getcdata_rev(nufid, arraysize, step, precision);
        if isempty(data)
            data = 0*u(:,:,kk);
        else
            data = reshape(data,Nc,Nkmax);
        end
        nuT(:,:,kk) = data;

        kk=kk+1;

    end

    % close files
    fclose all;
    clear sfid ufid Tfid etafid qfid

else % no data
    disp('no data files')
    return
end
    
%% convert EMPTY to nan

empties = u==EMPTY;
u(empties)=nan;
v(empties)=nan;
w(empties)=nan;
T(empties)=nan;
S(empties)=nan;
eta(eta==EMPTY)=nan;
q(empties)=nan;
nuT(empties)=nan;

% do same for AVG
if exist('AVG','var')
    AVG.empties = AVG.u==EMPTY;
    AVG.u(AVG.empties)=nan;
    AVG.v(AVG.empties)=nan;
    AVG.w(AVG.empties)=nan;
    AVG.T(AVG.empties)=nan;
    AVG.S(AVG.empties)=nan;
    AVG.eta(AVG.eta==EMPTY)=nan;
%     AVG.q(AVG.empties)=nan;
    AVG.nuT(AVG.empties)=nan;
%     AVG.u_avg(AVG.empties)=nan;
%     AVG.uw_var_avg(AVG.empties)=nan;
%     AVG.vw_var_avg(AVG.empties)=nan;
end


%% bottom cells
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

%% find densities
DATA = importdata('../state.c');
ttemp = strcmp('  return -prop->gamma*T;',DATA);
stemp = strcmp('  return prop->beta*s;',DATA);
tstemp = strcmp('  return prop->beta*s-prop->gamma*T;',DATA);

if max(ttemp)==1
    state='temp';
elseif max(stemp)==1
    state='sal';
elseif max(tstemp)==1
    state = 'full';
else
    state = 'sal';
    disp('caution - no equation of state match, assume salinity')
end
clear DATA stemp ttemp

if strcmp(state,'temp')
%     rho0 = rho0*(1 + beta.*nanmedian(reshape(S(:,:,1),[],1)));
    rho = rho0.*(1 -gamma.*T);
    if exist('AVG','var')
        AVG.rho = rho0.*(1 -gamma.*AVG.T);
    end
elseif strcmp(state,'sal')
%     rho0 = rho0*(1 - gamma.*nanmedian(reshape(T(:,:,1),[],1)));
    rho = rho0.*(1 + beta.*S);
    if exist('AVG','var')
        AVG.rho = rho0.*(1 + beta.*AVG.S);
    end
elseif strcmp(state,'full')
    rho = rho0.*(1+beta.*S-gamma.*T);
    if exist('AVG','var')
        AVG.rho = rho0.*(1+beta.*AVG.S-gamma.*AVG.T);
    end
end

% take rho background as avgerage distribution
rho_b = squeeze(nanmean(nanmean(rho(:,:,:),1),3))-rho0;
[rho_b,~,~]=meshgrid(rho_b,xv',t);
% take rhob as low pass filter of rho
% rho_b=nan+rho;
% for i=1:Nc
%     for j=1:Nkmax
%         if sum(~isnan(rho(i,j,:)))>5
%             [rho_b(i,j,:),~,~]=low_pass_filter(squeeze(rho(i,j,:))-rho0,dt,Tavg);
%         else
%             rho_b(i,j,:)=nan+  rho(i,j,:);
%         end
%     end
% end


% empties=isnan(u);
rho_b(empties)=nan;

rho_prime = rho - rho0 - rho_b;

% do same for AVG
if exist('AVG','var')
    AVG.rho0=rho0;
    AVG.rho_b = squeeze(nanmean(AVG.rho(:,:,1),1))-rho0;
    [AVG.rho_b,~,~]=meshgrid(AVG.rho_b,xv',AVG.t);
    AVG.rho_b(AVG.empties)=nan;
    AVG.rho_prime = AVG.rho - rho0 - AVG.rho_b;
end

%% pressures
% p = ph + q
% hydrostatic and dynamic pressures 
% ph = p0 + p_b(z) + p_prime(x,y,z,t)
% reference, background, and deviation

% reference hydrostatic pressure p0
i=1:Nc;
    for j=1:Nout_keep
        p0(i,:,j) = rho0*grav*(eta(i,:,j)-Z);
    end


% background pressure
p_b = zeros(Nc,Nkmax,Nout_keep);
p_b(:,1,:) = 0.5*rho_b(:,1,:).*grav.*(dz(1)+eta); % including the free surface in the pressure
%    here's a p_b calculation using cut cell dz values at the bottom
i=1:Nc;
    for k = 2:Nkmax
        p_b(i,k,:) = p_b(i,k-1,:)+0.5*rho_b(i,k-1,:).*grav.*dZ_cut(i,k-1)+...
                           0.5*rho_b(i,k,:).*grav.*dZ_cut(i,k);
    end
%set values below bottom to nan    
p_b = p_b +0*q;

% deviation pressure
p_prime = zeros(Nc,Nkmax,Nout_keep);
p_prime(:,1,:) = 0.5*rho_prime(:,1,:).*grav.*(dz(1)+eta); % including the free surface in the pressure
%    here's a p_prime calculation using cut cell dz values at the bottom
i=1:Nc;
    for k = 2:Nkmax
        p_prime(i,k,:) = p_prime(i,k-1,:)+0.5*rho_prime(i,k-1,:).*grav.*dZ_cut(i,k-1)+...
                           0.5*rho_prime(i,k,:).*grav.*dZ_cut(i,k);
    end
p_prime = p_prime +0*q;

% nonhydrostatic pressure
% add nonhydrostatic pressure to p_prime
p_hyd = p0 + p_b + p_prime;
% p_nhyd = q;
p_nhyd = q.*rho0;
p = p_hyd + p_nhyd;
% clear out any pressure in non-computational cells
p(empties==1) = nan;
p_hyd(empties==1) = nan;
p_nhyd(empties==1) = nan;
p_b(empties==1) = nan;
p_prime(empties==1) = nan;

%% save profile data

[PROF] = get_profile(datadir);

%% Brunt Vaisala Freq
% N^2 = -g/rho_0* d\rho/dz
% z vector is decreasing

N2 = nan+rho;
i=1:Nc;
    for j=1:Nout_keep
        N2(i,1:end-1,j) = -(grav./rho0) .* (rho(i,1:end-1,j) - rho(i,2:end,j)) ./ ...
            (0.5 .* (dZ(i,1:end-1) + dZ(i,2:end)));        
    end


%% Barotropic currents
for i=1:size(u,3)
    U(:,i) = nansum(u(:,:,i).*dZ_cut,2)./nansum(dZ_cut,2);
    V(:,i) = nansum(v(:,:,i).*dZ_cut,2)./nansum(dZ_cut,2);
%     pbar(:,i) = nansum(p_prime(:,:,i).*dZ_cut,2)./nansum(dZ_cut,2);
end
% baroclinic
for i=1:size(u,2)
    u_prime(:,i,:) = squeeze(u(:,i,:))-U;
    v_prime(:,i,:) = squeeze(v(:,i,:))-V;
%     p_prime(:,i,:) = squeeze(p_prime(:,i,:))-pbar;
end

% do same for AVG
if exist('AVG','var')
    for i=1:size(AVG.u,3)
        AVG.U(:,i) = nansum(AVG.u(:,:,i).*dZ_cut,2)./nansum(dZ_cut,2);
        AVG.V(:,i) = nansum(AVG.v(:,:,i).*dZ_cut,2)./nansum(dZ_cut,2);
    end
    % baroclinic
    for i=1:size(u,2)
        AVG.u_prime(:,i,:) = squeeze(AVG.u(:,i,:))-AVG.U;
        AVG.v_prime(:,i,:) = squeeze(AVG.v(:,i,:))-AVG.V;
    end
end



%% do time averaging, velocity decomposition
if exist('AVG','var')
    indxavg = t>(t(end)-Tavg);
    if length(AVG.t)>2
        indxAVG = AVG.t>(AVG.t(end)-Tavg+0.5*(AVG.t(2)));
    else
        indxAVG = length(AVG.t);
    end
    t_lowpass = t(end)-Tavg/2;
    % do bulk averaging
    % u = (Ubar + u_tildebar) + (Uprime+u_tildeprime);
    % mean barotropic, mean baroclinic, 
    % unsteady barotropic(tide), unsteady baroclinic (wave)

%     % U = Ubar+Uprime
%     Utilde = U - Ubar;
%     Vtilde = V - Vbar;
%     % u_tilde = u_tildebar+u_tildeprime
%     u_tildeprime = u_prime-u_barprime;
%     v_tildeprime = v_prime-v_barprime;

    Ubar = squeeze(nanmean(AVG.U(:,indxAVG),2));
    Vbar = squeeze(nanmean(AVG.V(:,indxAVG),2));
    
    u_barprime = nanmean(AVG.u_prime(:,:,indxAVG),3);
    v_barprime = nanmean(AVG.v_prime(:,:,indxAVG),3);
    
    for i=1:length(t)
       Utilde(:,i)=U(:,i)-Ubar;
       Vtilde(:,i)=V(:,i)-Vbar;
       
       u_tildeprime(:,:,i)=u_prime(:,:,i)-u_barprime;
       v_tildeprime(:,:,i)=v_prime(:,:,i)-v_barprime;
    end
    
    ubar = nanmean(AVG.u(:,:,indxAVG),3);
    vbar = nanmean(AVG.v(:,:,indxAVG),3);
    wbar = nanmean(AVG.w(:,:,indxAVG),3);    
    
    rhobar = nanmean(AVG.rho(:,:,indxAVG),3);
    etabar = nanmean(AVG.eta(:,:,indxAVG),3);
    nuTbar = nanmean(AVG.nuT(:,:,indxAVG),3);
    % don't need to save this
    rmfield(AVG,'nuT');

%     for i=1:Nc
%         Utilde(i,:)=U(i,:)-interp1(t_lowpass,Ubar(i,:),t');
%         Vtilde(i,:)=V(i,:)-interp1(t_lowpass,Vbar(i,:),t');
%         
% %         u_tildeprime(i,:,:) = 
%        for j=1:Nkmax
%            u_tildeprime(i,j,:)=squeeze(u_prime(i,j,:))'-interp1(t_lowpass,squeeze(u_barprime(i,j,:)),t');
%            v_tildeprime(i,j,:)=squeeze(v_prime(i,j,:))'-interp1(t_lowpass,squeeze(v_barprime(i,j,:)),t');
%         end
%      end

% do averaging at FFT
%      for i=1:Nc
%         [Ubar(i,:),~,~]=low_pass_filter(squeeze(U(i,:)),dt,Tavg);
%         [Vbar(i,:),~,~]=low_pass_filter(squeeze(V(i,:)),dt,Tavg);
%     end   
%     for i=1:Nc
%         for j=1:Nkmax
%             if sum(~isnan(u_prime(i,j,:)))>5
%                 % could do fancy FFT here
% %                 [u_barprime(i,j,:),~,~]=low_pass_filter(squeeze(u_prime(i,j,:)),dt,Tavg);
% %                 [v_barprime(i,j,:),~,~]=low_pass_filter(squeeze(v_prime(i,j,:)),dt,Tavg);
%             else
%               u_barprime(i,j,:)=nan+  u_prime(i,j,:);
%               v_barprime(i,j,:)=nan+  u_prime(i,j,:);
%             end
%               
%         end
%     end
    
    

    % u = (Ubar + u_tildebar) + (Uprime+u_tildeprime);
    % mean barotropic, mean baroclinic, 
    % unsteady barotropic(tide), unsteady baroclinic (wave)

%     % U = Ubar+Uprime
%     Utilde = U - Ubar;
%     Vtilde = V - Vbar;
%     % u_tilde = u_tildebar+u_tildeprime
%     u_tildeprime = u_prime-u_barprime;
%     v_tildeprime = v_prime-v_barprime;
% 
%     % interpolate down to save space
%     for i=1:Nc
%         temp1(i,:)=interp1(t,squeeze(Ubar(i,:)),t_lowpass);
%         temp2(i,:)=interp1(t,squeeze(Vbar(i,:)),t_lowpass);
%         for j=1:Nkmax
%             temp3(i,j,:)=interp1(t,squeeze(u_barprime(i,j,:)),t_lowpass);
%             temp4(i,j,:)=interp1(t,squeeze(v_barprime(i,j,:)),t_lowpass);
%         end
%     end
%     Ubar=temp1;
%     Vbar=temp2;
%     u_barprime=temp3;
%     v_barprime=temp4;
%     clear temp1 temp2 temp3 temp4
else % use instantaneous data for avg
    indxavg = t>(t(end)-Tavg);
    Ubar = nanmean(U(:,indxavg),2);
    Vbar = nanmean(V(:,indxavg),2);
    
    ubar = nanmean(u(:,:,indxavg),3);
    vbar = nanmean(v(:,:,indxavg),3);
    wbar = nanmean(w(:,:,indxavg),3);
    
    u_barprime = nanmean(u_prime(:,:,indxavg),3);
    v_barprime = nanmean(v_prime(:,:,indxavg),3);
    
     for i=1:length(t)
       Utilde(:,i)=U(:,i)-Ubar;
       Vtilde(:,i)=V(:,i)-Vbar;
       
       u_tildeprime(:,:,i)=u_prime(:,:,i)-u_barprime;
       v_tildeprime(:,:,i)=v_prime(:,:,i)-v_barprime;
    end
    
    rhobar = nanmean(rho(:,:,indxavg),3);
    etabar = nanmean(eta(:,indxavg),2);
    nuTbar = nanmean(nuT(:,:,indxavg),3);

end

%% Energy flux

Ek0 = 0.5*rho0*(U.^2+V.^2);
Ek_prime = 0.5*rho0*(u_prime.^2+v_prime.^2);
for i=1:Nkmax
    Ek0_prime(:,i,:) = rho0*(U.*squeeze(u_prime(:,i,:))...
                         +V.*squeeze(v_prime(:,i,:)));
end

% energy flux
for i=1:size(u,3)
    
    % barotropic energy flux, important terms are 2 and 3
    Fx_0(:,i) = U(:,i).*Ek0(:,i)+...
        rho0.*grav.*U(:,i).*depth.*squeeze(eta(:,1,i))+...
        U(:,i).*nansum(p_prime(:,:,i).*dZ,2)+...
        U(:,i).*nansum(p_nhyd(:,:,i).*dZ,2);
    
    Fy_0(:,i) = V(:,i).*Ek0(:,i)+...
        rho0.*grav.*V(:,i).*depth.*squeeze(eta(:,1,i))+...
        V(:,i).*nansum(p_prime(:,:,i).*dZ,2)+...
        V(:,i).*nansum(p_nhyd(:,:,i).*dZ,2);
    
    % baroclinic energy flux, important term is 3
%     Fx_prime(:,i) = nansum(u_tildeprime(:,:,i).*p_prime(:,:,i).*dZ,2); 
    Fx_prime(:,i) = nansum(u(:,:,i).*Ek_prime(:,:,i).*dZ,2)+...
        nansum(u(:,:,i).*Ek0_prime(:,:,i).*dZ,2)+...
        nansum(u_prime(:,:,i).*p_prime(:,:,i).*dZ,2)+...
        nansum(u_prime(:,:,i).*p_nhyd(:,:,i).*dZ,2);    
    
    Fy_prime(:,i) = nansum(v(:,:,i).*Ek_prime(:,:,i).*dZ,2)+...
        nansum(v(:,:,i).*Ek0_prime(:,:,i).*dZ,2)+...
        nansum(v_prime(:,:,i).*p_prime(:,:,i).*dZ,2)+...
        nansum(v_prime(:,:,i).*p_nhyd(:,:,i).*dZ,2);
end
%% simple averaging of last time step
indxavg = t>(t(end)-Tavg);

pbar = nanmean(p(:,:,indxavg),3);

Fx_0bar = nanmean(Fx_0(:,indxavg),2);
Fx_primebar = nanmean(Fx_prime(:,indxavg),2);

Fx_primebar_pos = Fx_prime(:,indxavg);
Fx_primebar_pos(Fx_primebar_pos<0)=nan;
Fx_primebar_pos = nanmean(Fx_primebar_pos,2);

Fx_primebar_neg = Fx_prime(:,indxavg);
Fx_primebar_neg(Fx_primebar_neg>0)=nan;
Fx_primebar_neg = nanmean(Fx_primebar_neg,2);

Fy_0bar = nanmean(Fy_0(:,indxavg),2);
Fy_primebar = nanmean(Fy_prime(:,indxavg),2);

%% get high pass variables

% indxavg2 = find(t>(t(end)-Tavg));
% thigh = t(indxavg);
% for i=1:length(indxavg2)
%     uhigh(:,:,i) = u(:,:,indxavg2(i)) - ubar;
%     vhigh(:,:,i) = v(:,:,indxavg2(i)) - vbar;
%     whigh(:,:,i) = w(:,:,indxavg2(i)) - wbar;
%     rhohigh(:,:,i) = rho(:,:,indxavg2(i)) - rhobar;
%     phigh(:,:,i) = p(:,:,indxavg2(i)) - pbar;
%     etahigh(:,:,i) = eta(:,indxavg2(i)) - etabar;
% end

%% do some stats on density
% [~,temp] = max(N2(:,:,1),[],2);
% indx_pycnocline(2) =round(nanmean(temp)); 
% D_pycnocline(2) = -z(round(nanmean(indx_pycnocline(2))));
% clear temp

rhotemp = nanmean(rho(:,:,1),1);

delta_rho = rhotemp(end)-rhotemp(1);

[~,indx_pycnocline(1)] = min(abs(quantile(rhotemp,0.3)-rhotemp));
[~,indx_pycnocline(2)] = min(abs(quantile(rhotemp,0.5)-rhotemp));
[~,indx_pycnocline(3)] = min(abs(quantile(rhotemp,0.7)-rhotemp));

D_pycnocline(1) = -z(indx_pycnocline(1));
D_pycnocline(2) = -z(indx_pycnocline(2));
D_pycnocline(3) = -z(indx_pycnocline(3));

%% track mid pycnocline 
rho_mid = nanmean(rho(:,indx_pycnocline(1),1))-rho0;
for i=1:size(rho,1)
    for j=1:size(rho,3)
        temp = squeeze(rho(i,:,j))-rho0;
        indx = ~isnan(temp) & ~isinf(temp); % make sure no nan or inf
        temp = temp(indx);
        tempz = z(indx);  
        % reorder to can interp
        [temp,ia,~]=unique(temp);
        tempz = tempz(ia);
              
        if length(temp)>2
            z_pyn(i,j) = interp1(temp,tempz,rho_mid)-z(indx_pycnocline(1));
        else
            z_pyn(i,j)=nan;
        end
    end
end
clear temp tempz ia indx

%% wavemaker results
if quad==0
clear indx r
Wave.x = x_wavemaker;
Wave.z = z';
Wave.t = t;
Wave.y = Ly/2;
Wave.z_pyc = z(indx_pycnocline(1));
% Wave.alphaw=AVG.alphaw;
if exist('AVG','var')
Wave.tavg=AVG.t;
end
for i=1:length(Wave.x) % find closest point
   r = sqrt((xv-Wave.x(i)).^2 + (yv-Wave.y).^2);
   [~, indx(i)] = min(r);
end
vmin = 0;
vmax = 0;
for i=1:length(Wave.x)
        Wave.u(i,:,:)=squeeze(u(indx(i),:,:));
        Wave.v(i,:,:)=squeeze(v(indx(i),:,:));
        Wave.w(i,:,:)=squeeze(w(indx(i),:,:));
        Wave.q(i,:,:)=squeeze(q(indx(i),:,:));
        
        
       
end

% find flux at b/c, midpoint, and just inside wavemakers
Wave.xF = [0 (0.2*sponge_distance) (sponge_distance)...
    (Lx/2) (Lx-sponge_distance) (Lx-0.2*sponge_distance) (Lx)];

for i=1:length(Wave.xF) % find closest point
   r = sqrt((xv-Wave.xF(i)).^2 + (yv-Wave.y).^2);
   [~, indx(i)] = min(r);
end
for i=1:length(Wave.xF)
    Wave.Fx_0(i)=Fx_0bar(indx(i));
    Wave.Fx_prime(i)=Fx_primebar(indx(i));
    Wave.Fx_prime_pos(i)=Fx_primebar_pos(indx(i));
    Wave.Fx_prime_neg(i)=Fx_primebar_neg(indx(i));
    Wave.uabs_max(i) = max(max(abs(u(indx(i),:,:))));
    Wave.uabs_pyc(i) = max(abs(u(indx(i),indx_pycnocline(1),:)));
    Wave.uabs_wave_pyc(i) = max(abs(u_tildeprime(indx(i),indx_pycnocline(1),:)));
    Wave.urms_wave_pyc(i) = sqrt(2*nanmean(u_tildeprime(indx(i),indx_pycnocline(1),:).^2));
end
end

%% box momentum balance

if quad
load SUNTANS_grid.mat
clear Quad

Quad.x = reshape(xv,Nx,Ny);
Quad.y = reshape(yv,Nx,Ny);
Quad.z = z;
Quad.dx = Quad.x(2,2)-Quad.x(1,1);
Quad.dy = Quad.y(2,2)-Quad.y(1,1);
Quad.dz = dz;
Quad.bottom_cells = reshape(bottom_cells,Nx,Ny);
Quad.depth = reshape(depth,Nx,Ny);
Quad.t = t;
Quad.rho0=rho0;
Quad.indxavg=indxavg;

Quad.etabar = reshape(etabar,Nx,Ny);
Quad.Ubar = reshape(Ubar,Nx,Ny);
Quad.Vbar = reshape(Vbar,Nx,Ny);

for i=1:length(z)
    for j=1:length(t)
        Quad.u(:,:,i,j) = reshape(u(:,i,j),Nx,Ny);
        Quad.v(:,:,i,j) = reshape(v(:,i,j),Nx,Ny); 
        Quad.w(:,:,i,j) = reshape(w(:,i,j),Nx,Ny); 
        Quad.rho(:,:,i,j) = reshape(rho(:,i,j),Nx,Ny);
        Quad.p_hyd(:,:,i,j) = reshape(p_hyd(:,i,j),Nx,Ny);
        Quad.p_nhyd(:,:,i,j) = reshape(p_nhyd(:,i,j),Nx,Ny);
        
        Quad.p_b(:,:,i,j) = reshape(p_b(:,i,j),Nx,Ny); 
        Quad.p0(:,:,i,j) = reshape(p0(:,i,j),Nx,Ny); 
        Quad.p_prime(:,:,i,j) = reshape(p_prime(:,i,j),Nx,Ny); 
        
        Quad.nuT(:,:,i,j) = reshape(nuT(:,i,j),Nx,Ny);
        
    end
    Quad.ubar(:,:,i) = reshape(ubar(:,i),Nx,Ny);
    Quad.vbar(:,:,i) = reshape(vbar(:,i),Nx,Ny);
    Quad.wbar(:,:,i) = reshape(wbar(:,i),Nx,Ny);
    Quad.pbar(:,:,i) = reshape(pbar(:,i),Nx,Ny);
    
    Quad.nuTbar(:,:,i) = reshape(nuTbar(:,i),Nx,Ny);

end
for j=1:length(t)
   Quad.U(:,:,j) = reshape(U(:,j),Nx,Ny);
   Quad.V(:,:,j) = reshape(V(:,j),Nx,Ny);
   Quad.eta(:,:,j) = reshape(eta(:,j),Nx,Ny);
end

for i=1:length(t) % spatial average
Quad.U_avg(i) = sum(reshape(Quad.U(:,:,i).*Quad.depth,[],1))./...
    sum(reshape(Quad.depth,[],1));
Quad.U_std(i) = nanstd(reshape(Quad.U(:,:,i),[],1));

Quad.V_avg(i) = sum(reshape(Quad.V(:,:,i).*Quad.depth,[],1))./...
    sum(reshape(Quad.depth,[],1));
Quad.V_std(i) = nanstd(reshape(Quad.V(:,:,i),[],1));

end

% dynamic pressure = non-hydrostatic + free surface
for i=1:Nx
    for j=1:Ny
        for k=1:Nout_keep
            Quad.p_dyn(i,j,:,k) = Quad.p_nhyd(i,j,:,k) + rho0*grav*Quad.eta(i,j,k);
        end
    end
end
%% avg vorticity

Quad.vorty_bar = nan+Quad.ubar;

for i=2:Nx-1
    for j=2:Nkmax-1
        Quad.vorty_bar(i,:,j) = -((Quad.ubar(i,:,j-1)-Quad.ubar(i,:,j+1))/(Quad.z(j-1)-Quad.z(j+1))...
            -(Quad.wbar(i+1,:,j)-Quad.wbar(i-1,:,j))./(2*dx));
    end
end

%% compute spatial profile



%% periodic bc calculation of Cd
f1 = dir('../data/KurtS.dat');
f2 = importdata('../sources.c');
iforce = strcmp('REAL FORCE;',f2);
if ~isempty(f1)
    % units are [m/s^2]
    Saccel = importdata([f1.folder '/' f1.name]);
    [Saccelt,ia,ic] = unique(Saccel(:,2)');
    DRAG.Saccel = Saccel(ia,1);
    DRAG.Saccel = DRAG.Saccel(1:istep:end);
elseif max(iforce)==1
f = [];    j=1; 
    while isempty(f) && j<length(f2)
        [f] = sscanf(f2{j},'REAL ULm = %f, Cd=%f, H=%f, Lx=%f, TL = %f;');
        j=j+1;
    end
    if isempty(f)
        disp('could not read sources.c')
        Ulm=nan;
        Cd=nan;
        H=nan;
    end
    Ulm=f(1);
    Cd=f(2);
    H=f(3);
    DRAG.Saccel = Cd*Ulm.^2/H+0*t;
    
else
    DRAG.Saccel=nan+t;
end
clear DATA iforce f j Ulm Cd H f1 f2 Saccelt ia ic



%% bulk momentum balance for non periodic B/C
% fit pg to first 80% of slope

xdata = xv;
ydata = squeeze(eta(:,:,end));
indx = xdata< 0.8*L & ~isnan(ydata);
xdata = xdata(indx);
ydata = ydata(indx);

pfit = polyfit(xdata,ydata,1);

Quad.xfit = linspace(min(xv),max(xv));
Quad.yfit = polyval(pfit,Quad.xfit);


% plot(xv,squeeze(eta(:,:,end)),'.',...
%     xfit,yfit,'-')

Quad.PG_avg = pfit(1);
Quad.uabsu_h = (Quad.U(:,:,end).*sqrt(Quad.U(:,:,end).^2+Quad.V(:,:,end).^2))./...
        (Quad.depth+Quad.eta(:,:,end));
Quad.uabsu_h = nanmean(reshape(Quad.uabsu_h,[],1));    

Quad.Cd_avg = -grav*Quad.PG_avg./Quad.uabsu_h;



%% compute spatial variation of drag coeff

% get bulk <Cd>
DRAG.Ub = nanmean(Quad.U_avg(indxavg));


DRAG.rho = rho0;
DRAG.H = nanmean(depth+etabar);

% save total force [N] = Saccel * \rho * Vol
DRAG.ForceB = DRAG.Saccel*rho0*DRAG.H*L*W;

indxavg = t>(t(end)-Tavg);

for i=1:length(DRAG.Saccel)
    DRAG.Cd_bulk(:,:,i) = 0*Quad.depth+DRAG.H.*DRAG.Saccel(i)./(DRAG.Ub*abs(DRAG.Ub));
%     DRAG.Cd_bulk(:,:,i) = Quad.depth.*DRAG.Saccel(i)./(Quad.U(:,:,i).*abs(Quad.U(:,:,i)));
end


DRAG.Cd_bulk_bar = nanmean(DRAG.Cd_bulk(:,:,Quad.indxavg),3);
DRAG.Cd_bulk_avg = squeeze(nanmean(nanmean(DRAG.Cd_bulk,1),2));
DRAG.Cd_bulk_avg_bar = nanmean(DRAG.Cd_bulk_avg(Quad.indxavg));

% get bottom dz, u, v
for i=1:Nx
    for j=1:Ny
        Quad.dzB(i,j) = Quad.dz(Quad.bottom_cells(i,j));
        
        Quad.ubot(i,j,:) = squeeze(Quad.u(i,j,Quad.bottom_cells(i,j),:));
        Quad.vbot(i,j,:) = squeeze(Quad.v(i,j,Quad.bottom_cells(i,j),:));
        
        Quad.pbot(i,j,:) = squeeze(Quad.p_dyn(i,j,Quad.bottom_cells(i,j),:));
        
    end
end

% get slope dhdx, get depth on edges and subtract
Quad.zbot = -Quad.depth;
% Quad.dhdx=nan+Quad.zbot;
% for i=1:Nx-1
%      for j=1:Ny
%          Quad.dhdx(i,j) = (Quad.zbot(i+1,j)-Quad.zbot(i,j))/dx;
%      end
% end


% get drag coefficient
% if nu==0 % no viscous effects
%     DRAG.CdB=0*Quad.depth;
if CdB==0 % specified z0
    DRAG.z0B=z0B;
    DRAG.CdB = 0.42^2*(log((0.5*Quad.dzB) ./ (z0B))).^-2;
else % specifified Cd
    DRAG.CdB = CdB + 0*Quad.depth;
end

% bottom shear stress
for i=1:Nout_keep
    DRAG.taux(:,:,i) = DRAG.rho.*DRAG.CdB.*Quad.ubot(:,:,i).*...
        sqrt(Quad.ubot(:,:,i).^2+Quad.vbot(:,:,i).^2);
    DRAG.tauy(:,:,i) = DRAG.rho.*DRAG.CdB.*Quad.vbot(:,:,i).*...
        sqrt(Quad.ubot(:,:,i).^2+Quad.vbot(:,:,i).^2);
end

DRAG.Cd_taux = DRAG.taux./(DRAG.rho*DRAG.Ub*abs(DRAG.Ub));
DRAG.Cd_taux_bar = nanmean(DRAG.Cd_taux(:,:,indxavg),3);
DRAG.Cd_taux_avg = squeeze(nanmean(nanmean(DRAG.Cd_taux,1),2));
DRAG.Cd_taux_avg_bar = nanmean(DRAG.Cd_taux_avg(indxavg));

%% pressure drag

% get heights of cells
Quad.zb = Quad.z-dz/2; % bottom face
Quad.zt = Quad.z+dz/2; % top face
for i=1:Nx
    for j=1:Ny
        % cells below depth = 0
        Quad.h(i,j,:) = nan;
        % cells fully open
        indx = Quad.zb>=-Quad.depth(i,j);
        Quad.h(i,j,indx)= dz(indx);
        % cells partially closed
        indx = Quad.zb<-Quad.depth(i,j) & Quad.zt>-Quad.depth(i,j);
        Quad.h(i,j,indx)= Quad.zt(indx)+Quad.depth(i,j);
    end
end

% now compute the differential height of each face
Quad.dhface = nan+Quad.h;
for i=1:Nx-1 % first all partial cells & cells with no face height
    j=1:Ny;
    k=1:Nkmax;
    Quad.dhface(i,j,k) =Quad.h(i+1,j,k)-Quad.h(i,j,k);
end

% now compute face heights for fully blocked + is on R, - is on left
for i=1:Nx-1
    for j=1:Ny
        for k=1:Nkmax
            if isnan(Quad.h(i,j,k)) && isnan(Quad.h(i+1,j,k)) % no face,fully underground
                Quad.dhface(i,j,k)=nan;
            elseif ~isnan(Quad.h(i,j,k)) && isnan(Quad.h(i+1,j,k)) % right face
                Quad.dhface(i,j,k)=+Quad.h(i,j,k);
            elseif isnan(Quad.h(i,j,k)) && ~isnan(Quad.h(i+1,j,k)) % left face
                Quad.dhface(i,j,k) = -Quad.h(i+1,j,k);
            
            end
        end
    end
end
% pcolorjw(flipud(squeeze(Quad.dhface(:,4,:,end))')), colorbar


% take pressure gradient on edges

% find pressure on face of cells,
% indexing is face i is Right side of cell i

Quad.pface = nan+Quad.p_dyn; % sum of q and rho*g*eta
for i=1:Nx-1
   Quad.pface(i,:,:,:)=0.5*(Quad.p_dyn(i,:,:,:)+Quad.p_dyn(i+1,:,:,:)); 
end
% last cell due to periodicity
Quad.pface(end,:,:,:)=0.5*(Quad.p_dyn(1,:,:,:)+Quad.p_dyn(end,:,:,:));

% pface is good for all fully open cells
% pface at this point is nan for any face touching full cell block
% need to interp pface from dp/dx from left and right
for i=2:Nx-2
    for j=1:Ny
        for k=1:Nkmax
            if Quad.dhface(i,j,k)>0 % try interp from left pf(i)=pcell(i)+(pcell(i)-pcell(i-1))/(0.5dx)
%                 Quad.pface(i,j,k,:)=Quad.p_dyn(i,j,k,:);
                  Quad.pface(i,j,k,:)=Quad.p_dyn(i,j,k,:)+...
                    (Quad.p_dyn(i,j,k,:)-Quad.p_dyn(i-1,j,k,:))/(0.5*Quad.dx);
            end
            if Quad.dhface(i,j,k)<0% try interp from right pf(i)=pcell(i)-(pcell(i+1)-pcell(i))/(0.5dx)
%                 Quad.pface(i,j,k,:)=Quad.p_dyn(i+1,j,k,:);
                  Quad.pface(i,j,k,:)=Quad.p_dyn(i+1,j,k,:)-...
                    (Quad.p_dyn(i+2,j,k,:)-Quad.p_dyn(i+1,j,k,:))/(0.5*Quad.dx);
            end
        end
    end
end

% pcolorjw(flipud(squeeze(Quad.pface(:,4,:,end))')), colorbar


% now compute volume averaged acceleration from pressure forces
% F = m * a;
% a = F ./ m = F ./ (\rho * Vol_cell)
% acc_px = sum_z(pface * dhf) / (\rho * sum_z(dhface) * dx)
for i=1:length(Quad.t)

    DRAG.acc_px(:,:,i) = nansum(Quad.pface(:,:,:,i).*Quad.dhface,3) ./...
        (DRAG.rho*nansum(abs(Quad.h),3)*Quad.dx);
end

% for i=1:Nout_keep
%     DRAG.Cd_px(:,:,i) = Quad.depth.*DRAG.acc_px(:,:,i)./(Quad.U(:,:,i).*abs(Quad.U(:,:,i)));
% end
DRAG.Cd_px = DRAG.acc_px.*DRAG.H/(DRAG.Ub.^2);

DRAG.Cd_px_bar = nanmean(DRAG.Cd_px(:,:,Quad.indxavg),3);
DRAG.Cd_px_avg = squeeze(nanmean(nanmean(DRAG.Cd_px,1),2));
DRAG.Cd_px_avg_bar = nanmean(DRAG.Cd_px_avg(Quad.indxavg));

%% convert Cd to z0 equivalent

DRAG.z0_bulk = DRAG.H ./ exp(0.41/sqrt(DRAG.Cd_bulk_avg_bar)-(0.2-1));
DRAG.z0_px = DRAG.H ./ exp(0.41/sqrt(DRAG.Cd_px_avg_bar)-(0.2-1));
DRAG.z0_taux = DRAG.H ./ exp(0.41/sqrt(DRAG.Cd_taux_avg_bar)-(0.2-1));


%% Reynolds stresses
up = Quad.u(:,:,:,indxavg);
vp = Quad.v(:,:,:,indxavg);
wp = Quad.w(:,:,:,indxavg);
for i=1:size(up,4)
   up(:,:,:,i) = up(:,:,:,i)-nanmean(up,4);
   vp(:,:,:,i) = vp(:,:,:,i)-nanmean(vp,4);
   wp(:,:,:,i) = wp(:,:,:,i)-nanmean(wp,4);
   
end
Quad.upup_bar = nanmean(up.*up,4);
Quad.upvp_bar = nanmean(up.*vp,4);
Quad.upwp_bar = nanmean(up.*wp,4);
Quad.vpvp_bar = nanmean(vp.*vp,4);
Quad.vpwp_bar = nanmean(vp.*wp,4);
Quad.wpwp_bar = nanmean(wp.*wp,4);
clear up vp wp

%% find Reynolds stresses on bottom
for i=1:size(Quad.x,1)
    for j=1:size(Quad.x,2)
        Quad.upwp_bot(i,j) = Quad.upwp_bar(i,j,Quad.bottom_cells(i,j));
        Quad.upup_bot(i,j) = Quad.upup_bar(i,j,Quad.bottom_cells(i,j));
    end
end
DRAG.Cd_Re_bar = -Quad.upwp_bot ./ DRAG.Ub^2;
DRAG.Cd_Re_avg_bar = nanmean(reshape(DRAG.Cd_Re_bar,[],1));

%% 2D momentum balance
grav=9.81;

% all terms on LHS
% US + NL1 + NL2 + PG - Fs + BT + Fp + HD + Re + Err= 0;

% compute velocities on faces
% center cells
Quad.Uface = 0.5*(Quad.U(2:end,:,end)+Quad.U(1:end-1,:,end));
% periodic contraint, add in end faces (which are the same)
utemp = 0.5*(Quad.U(1,:,end)+Quad.U(end,:,end));
Quad.Uface = [utemp; Quad.Uface; utemp];

% center cells
Quad.Vface = 0.5*(Quad.V(:,2:end,end)+Quad.V(:,1:end-1,end));
% boundary conditions set v=0 at walls
utemp = 0*Quad.Vface(:,1);
Quad.Vface = [utemp, Quad.Vface, utemp];

% compute terms at cell centers
% Unsteady
Quad.USx = (Quad.U(:,:,end) - Quad.U(:,:,end-1))./...
            (Quad.t(end)-Quad.t(end-1));
Quad.USx_bar = nanmean(nanmean(abs(Quad.USx)));

% nonlinearity 1
Quad.NL1 = Quad.U(:,:,end).*(Quad.Uface(2:end,:) - Quad.Uface(1:end-1,:))./...
       (Quad.dx);
% Quad.NL1(2:end,:) = (Quad.U(2:end,:,end).^2 - Quad.U(1:end-1,:,end).^2)./...
%        (Quad.x(2:end,:) - Quad.x(1:end-1,:));
Quad.NL1_bar = nanmean(nanmean(abs(Quad.NL1)));

% nonlinearity 2, use CD for U derivative
Quad.NL2 = nan+Quad.NL1;
Quad.NL2(:,2:end-1) = Quad.V(:,2:end-1,end).*(Quad.U(:,3:end,end) - Quad.U(:,1:end-2,end))./...
       (2*Quad.dy);
Quad.NL2_bar = nanmean(nanmean(abs(Quad.NL2)));

% Pressure gradient
etemp = Quad.eta(:,:,end);
% add extra row of x for periodicity
etemp = [etemp; etemp(1,:)];
Quad.PGx = grav*(etemp(2:end,:) - etemp(1:end-1,:,end))./...
       (Quad.dx);
Quad.PGx_bar = nanmean(nanmean(abs(Quad.PGx)));

% Driving pressure force S
Quad.Fs = -DRAG.Saccel(end) + 0*Quad.USx;
Quad.Fs_bar = -nanmean(nanmean(abs(Quad.Fs)));

% bottom stress
Quad.BT = DRAG.taux(:,:,end) ./ (Quad.depth+Quad.eta(:,:,end))./rho0;
Quad.BT_bar = nanmean(nanmean(abs(Quad.BT)));

% Reynolds shear stress
Quad.ReS = Quad.upwp_bot ./ (Quad.depth+Quad.eta(:,:,end));
Quad.ReS_bar = nanmean(nanmean(abs(Quad.ReS)));

% bottom pressure force
Quad.Fpx = DRAG.acc_px(:,:,end);
Quad.Fpx_bar = nanmean(nanmean(abs(Quad.Fpx)));

% horizontal diffusion, use central difference to get d2U/dx2
% nutemp = nanmean(Quad.nuT(:,:,:,end),3);
Quad.HD = nan+Quad.USx;
Quad.HD(2:end-1,:) = -nu_H.*(Quad.U(3:end,:,end) +...
    2*Quad.U(2:end-1,:,end) + Quad.U(1:end-2,:,end))/Quad.dx^2;
Quad.HD_bar = nanmean(nanmean(abs(Quad.HD)));

% residual error
Quad.Err = Quad.USx + Quad.NL1 + Quad.NL2 + Quad.PGx + Quad.Fs + Quad.BT + Quad.Fpx + Quad.HD + Quad.ReS;
Quad.Err_bar = nanmean(nanmean(abs(Quad.Err)));

% Quad.Cd = -Quad.PGx .*(Quad.depth+Quad.eta(:,:,end))./...
%     (Quad.U(:,:,end).*sqrt(Quad.U(:,:,end).^2+Quad.V(:,:,end).^2));
% 
% Quad.Cd_bar = nanmedian(reshape(Quad.Cd,[],1));
clear etemp utemp nutemp

%% average profile
DRAG.PROF.z = Quad.z'-min(Quad.z);
for i=1:length(Quad.z)
    temp = reshape(Quad.u(:,:,i,Quad.indxavg),[],1);
    DRAG.PROF.ubar(i) = nanmean(temp);
    DRAG.PROF.ustd(i) = nanstd(temp);
end
DRAG.PROF.ztop = -min(reshape(Quad.depth,[],1))-min(Quad.z);
DRAG.PROF.zmid = -nanmean(reshape(Quad.depth,[],1))-min(Quad.z);
DRAG.PROF.Ub = DRAG.Ub;



%% stats on bathymetry
clear BATH 

BATH.x = Quad.x;
BATH.y = Quad.y;
BATH.z = Quad.z;
BATH.Z = -Quad.depth;
BATH.D0 = nanmean(reshape(Quad.depth,[],1));

BATH.dx = BATH.x(2,2)-BATH.x(1,1);
BATH.dy = BATH.y(2,2)-BATH.y(1,1);

BATH.dzdx = nan+BATH.x;
BATH.dzdx(2:end,:) = (BATH.Z(2:end,:) - BATH.Z(1:end-1,:))/BATH.dx;

BATH.dzdy = nan+BATH.y;
BATH.dzdy(:,2:end) = (BATH.Z(:,2:end) - BATH.Z(:,1:end-1))/BATH.dy;

BATH.Ub = DRAG.Ub;
BATH.nu = 1e-6;

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

% get power spectral density
BATH.Szx = BATH.Fzx.*conj(BATH.Fzx)./(BATH.dfx);
BATH.Sdzdx = BATH.Fdzdx.*conj(BATH.Fdzdx)./(BATH.dfx);

BATH.Szy = BATH.Fzy.*conj(BATH.Fzy)./(BATH.dfy);
BATH.Sdzdy = BATH.Fdzdy.*conj(BATH.Fdzdy)./(BATH.dfy);

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
BATH.Re_h = BATH.Ub*BATH.hrms/BATH.nu;
BATH.h_lambdax = sqrt(2)/pi*BATH.dzdx_rms;
if size(BATH.x,2)>3
    BATH.h_lambday = sqrt(2)/pi*BATH.dzdy_rms;
else
    BATH.h_lambday=0;
end

%% save Quad format file
close all

save([savedir 'SUNTANS_quad_results.mat'],'Quad','DRAG','BATH','-v7.3')

else % keep in unstructured format 
% clean up and save
disp('saving MAT file')

clear ans arraysize bot c cols data F i j k is ntout precision step tend teval ...
    Quad BATH DRAG Ek0_prime Ek_prime Uface empties p0 tstemp

save([savedir 'SUNTANS_results.mat'],'-v7.3')

%% save only last time step and save small file
clear AVG indxavg
eta = eta(:,:,end);
Ek0 = Ek0(:,end);
% Fx_0 = Fx_0(:,end);
% Fx_prime = Fx_prime(:,end);
% Fy_0 = Fy_0(:,end);
% Fy_prime = Fy_prime(:,end);
N2 = N2(:,:,end);
nuT = nuT(:,:,end);
p = p(:,:,end);
p_b = p_b(:,:,end);
p_hyd = p_hyd(:,:,end);
p_nhyd = p_nhyd(:,:,end);
p_prime = p_prime(:,:,end);
q = q(:,:,end);
rho = rho(:,:,end);
rho_b = rho_b(:,:,end);
rho_prime = rho_prime(:,:,end);
S = S(:,:,end);
T = T(:,:,end);
t_ref = t_ref(end);
u = u(:,:,end);
U = U(:,end);
u_prime = u_prime(:,:,end);
u_tildeprime = u_tildeprime(:,:,end);
Utilde = Utilde(:,end);
v = v(:,:,end);
V = V(:,end);
v_prime = v_prime(:,:,end);
v_tildeprime = v_tildeprime(:,:,end);
Vtilde = Vtilde(:,end);
w = w(:,:,end);
z_pyn = z_pyn(:,end);

disp('saving end MAT file')

save([savedir 'SUNTANS_results_end.mat'],'-v7.3')
end   


%% stuff from Eric
% for n=1:nout
%     udata=fread(ufid,Nc*Nkmax*3,'float64');
%     ud = reshape(udata,Nc,Nkmax,3);
%     u = squeeze(ud(:,:,1));
%     w = squeeze(ud(:,:,3));
%     Tdata=fread(Tfid,Nc*Nkmax,'float64');
%     T = reshape(Tdata,Nc,Nkmax);
%     eta = fread(eta_fid,Nc,'float64');
%     
%     % for the first nout, define bottom cells, fix dZ in bottom cells
%     % and the background temperature profile, T_bar
%     if(n==1)
%         empties = zeros(size(u));
%         empties(u==EMPTY)=1;
%         u(u==EMPTY)=nan;
%         w(w==EMPTY)=nan;
%         T(T==EMPTY)=nan;
%         bottom_cells = zeros(Nc,1);
%             for i = 1:Nc
%                 if isnan(T(i,Nkmax))
%                     bottom_cells(i) = find(isnan(T(i,:)),1)-1;
%                 else
%                     bottom_cells(i) = Nkmax;
%                 end
%                 bot = bottom_cells(i);
%                 dZ(i,bot) = depth(i)+(z(bot)+dz(bot)/2);
%             end  
%         % fix any empty cells now to 0
%         u(isnan(u))=0;
%         w(isnan(w))=0;
%         T(isnan(T))=0;
%         
%         % save this w for animation
%         wanimate(n,:,:) = w;
% 
%         %compute initial temperature profile
%         T0 = 20 + (N^2/(gamma*grav)).*z'; %computed using N, as was done in SUNTANS
%         Tdiff = T(1,:)-T0; %compare initial data T to this theoretical setting (should be spot on)
%         %     figure();
%         %     plot(Tdiff,z)
%         teval = find(depth == max(depth),1); % make sure the T_bar taken from the data is from the deepest place!
%         T_bar=ones(Nc,1)*T(teval,:);
%         %T_bar = ones(Nc,1)*T0;
%     end
%     
%     
%     % if n ~=1, still have to fix any empty cells to 0's
%     u(u==EMPTY)=0;
%     w(w==EMPTY)=0;
%     T(T==EMPTY)=0;
%     wanimate(n,:,:) = w;
%     
%     % density
%     rho = rho0.*(-gamma.*T);
%     rho_prime = -gamma*rho0*(T-T_bar);
%     
%     % hydrostatic pressure (with/without use of free surface
%     % deflections and/or cut cell bottm dz's)
%     p_prime = zeros(Nc,Nkmax);
%     p_prime(:,1) = 0.5*rho_prime(:,1)*grav.*(dz(1)+2*eta.*1e8); % including the free surface in the pressure
% %     p_prime(:,1) = 0.5*rho_prime(:,1)*grav.*(dz(1));
% %     for k=2:Nkmax
% %         p_prime(:,k) = p_prime(:,k-1)+0.5*rho_prime(:,k-1)*grav*dz(k-1)+...
% %                                       0.5*rho_prime(:,k)*grav*dz(k);
% %     end
% %    here's a p_prime calculation using cut cell dz values at the bottom
%     for i=1:Nc
%         for k = 2:bottom_cells(i)
%         p_prime(i,k) = p_prime(i,k-1)+0.5*rho_prime(i,k-1)*grav*dZ(i,k-1)+...
%                                       0.5*rho_prime(i,k)*grav*dZ(i,k);
%         end
%     end
%     
%     % nonhydrostatic pressure
%     q = reshape(fread(qfid,Nc*Nkmax,'float64'),Nc,Nkmax);
%     % add nonhydrostatic pressure to p_prime
%     p_prime = p_prime+q.*rho0;
%     % clear out any pressure in non-computational cells
%     p_prime(empties==1) = 0;
%     
%     % fill animation arrays
%     panimate(n,:,:) = p_prime;
%     tanimate(n,:,:) = T;
% 
%     % Barotropic currents
%     U = sum(u.*dZ,2)./sum(dZ,2)*ones(1,Nkmax);
%      
%     % compute p_prime at the floor (rather than cell center in bottom cell) 
%     % and calculate form drag related quantities 
%     for i = 1:Nc
%         kb = bottom_cells(i);
% %         dZ(i,kb) = dZ(i,kb) - h(i) + (D+z(kb)-0.5*dz(kb));
%         
%         % p_prime_bottom = p_prime(i,kb-1) + rho_prime(i,kb)*dZ(i,kb); 
%         p_prime_bottom = p_prime(i,kb) + rho_prime(i,kb)*0.5*dZ(i,kb); 
%         
%         center_force_per_area = p_prime(i,kb)*dhdx(i);
%         force_per_area = p_prime_bottom*dhdx(i);
%         power_per_area = U(i,kb)*force_per_area;
%         force_per_y(n) = force_per_y(n) + force_per_area * dx;
%         power_per_y(n) = power_per_y(n) + power_per_area * dx;
%         cell_center_power(n) = cell_center_power(n) + U(i,kb)*center_force_per_area*dx*dy;
%     end
% 
%     pw = p_prime.*w.*dx*dy;
%     pu = p_prime.*(u-U).*dZ.*dy;
%     Fv = sum(pw,1);
%     Fh = sum(pu,2);
% 
%     % Box flux
%     F_l(n) = sum(pu(i1,k_flux:Nkmax));
%     F_r(n) = sum(pu(i2,k_flux:Nkmax));
%     F_t(n) = sum(pw(i1:i2,k_flux));    
%     
%     %  add to pw_on_hill
%     pw_on_hill(:,n) = squeeze(p_prime(floor(Nc/2),:).*w(floor(Nc/2),:));
% 
%     % energy calculations
%     Ek(n) = sum(sum((u.*u + w.*w).*rho.*dV));
%     Ep(n) = sum(sum(grav.*rho.*(D+Z).*dV));
%     Etot(n) = Ek(n)+Ep(n);
%     Ekprime(n) = sum(sum(((u-U).*(u-U) + w.*w).*rho.*dV));
%     Epprime(n) = sum(sum(grav.*rho_prime.*(D+Z).*dV));
%     
%     
%     % plotting p' and fluxes at three different times in run
%     if (makeplots==1 && (n== floor(nout/3) || n== floor(nout/2) || n==nout))   
%         figure()
%         subplot(3,2,[1,2])
%         contourf(X,Z,p_prime,'LineStyle','none')
%         colorbar
%         title(['perturbation pressure (kg/m/s^2) at t=',num2str(n*ntout*dt/tday,3),' days'])
%         xlabel('easting (m)')
%         ylabel('depth (m)')
%         hold on
%         plot((xv(i1)-xc).*ones(size(k_flux:Nkmax)),z(k_flux:Nkmax),'r',...
%             xv(i1:i2)-xc,z(k_flux).*ones(size(i1:i2)),'r',...
%             (xv(i2)-xc).*ones(size(k_flux:Nkmax)),z(k_flux:Nkmax),'r')
%         hold off
%         
%         subplot(3,2,[3,4])
%         plot(xv(i1:i2)-xc,pw(i1:i2,k_flux)./(pw_Gill*dx*dy))
%         title(['p*w/pw_{Gill} along top of box, z=',num2str(D+z(k_flux)),'m'])
%         xlabel('easting (m)')
%         ylabel('p*w/pw_{Gill}')
%         
%         subplot(3,2,5)
%         plot(pu(i1,k_flux:Nkmax)./(dZ(k_flux:Nkmax).*dy*pw_Gill),D+z(k_flux:Nkmax))
%         title(['p*u along left of box, x=',num2str(xv(i1)-xc),'m'])
%         xlabel('p*u/pw_{Gill}')
%         ylabel('depth (m)')
%         
%         subplot(3,2,6)
%         plot(pu(i2,k_flux:Nkmax)./(dZ(k_flux:Nkmax).*dy*pw_Gill),D+z(k_flux:Nkmax))
%         title(['p*u along right of box, x=',num2str(xv(i2)-xc),'m'])
%         xlabel('p*u/pw_{Gill}')
%         ylabel('depth (m)')
%         
%     end
%     
% end
% 
% force_from_hill = force_per_y.*dy;
% power_from_hill = power_per_y.*dy;
% 
% %% Saving
% save([savedir,savename])
% 
% %% Plotting energy flux vs time
% if makeplots==1
%     figure()  
%     subplot(2,1,1)
%     plot(t/tday,F_r/E_Gill,'k:',...
%          t/tday,F_l/E_Gill,'k-.',...
%          t/tday,F_t/E_Gill,'k--',...
%          t/tday,(F_r-F_l+F_t)/E_Gill,'k-');
%     title('Fluxes out of a box around the hill');
%     legend(sprintf('Right x=+%.1f km',(xv(i1)-xc)/1000),...
%            sprintf('Left x=%.1f km',(xv(i2)-xc)/1000),...
%            sprintf('Top z=%.1f m',h_flux),...
%            'Net (Top + Right - Left)',...
%            'location','east');
%     xlabel('Time (d)');
%     ylabel('Flux/F_{Gill} @ z_{flux}');
% 
%     subplot(2,1,2)
%     plot(t/tday,power_from_hill/E_Gill,'k.-',...
%         t/tday,cell_center_power/E_Gill,'k:',...
%         t/tday,ones(size(t)),'k--');
%     title('Flux via form drag');
%     xlabel('Time (d)');
%     ylabel('Flux/F_{Gill} @ z_{flux}');
%     legend('p at cell bottom','p at cell center','location','east')
% end
% 
% %% Cleaning up
% fclose all;
% 
