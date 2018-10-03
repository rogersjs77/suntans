clear
% load and save, surface results for movie

load('./SUNTANS_results.mat',...
    'xv','yv','z','x_sponge','t','Nkmax',...
    'u_tildeprime','v_tildeprime','ubar','vbar','eta',...
    'T','datadir','indx_pycnocline','z_pyn','rho','rho0');

GRID = load(['./SUNTANS_grid.mat']);

% %% HACK - recompute z_pyn
% rho_mid = nanmean(rho(:,indx_pycnocline(1),1))-rho0;
% for i=1:size(rho,1)
%     for j=1:size(rho,3)
%         temp = squeeze(rho(i,:,j))-rho0;
%         indx = ~isnan(temp) & ~isinf(temp); % make sure no nan or inf
%         temp = temp(indx);
%         tempz = z(indx);  
%         % reorder to can interp
%         [temp,ia,~]=unique(temp);
%         tempz = tempz(ia);
%               
%         if length(temp)>2
%             z_pyn(i,j) = interp1(temp,tempz,rho_mid)-z(indx_pycnocline(2));
%         else
%             z_pyn(i,j)=nan;
%         end
%     end
% end
% clear temp tempz ia indx rho
%%


eta = reshape(eta,GRID.Nx,GRID.Ny,[]);
z_pyn = reshape(z_pyn,GRID.Nx,GRID.Ny,[]);

u_tildeprime_depth = 100*reshape(u_tildeprime(:,indx_pycnocline(1),:),GRID.Nx,GRID.Ny,[]);
v_tildeprime_depth = 100*reshape(v_tildeprime(:,indx_pycnocline(1),:),GRID.Nx,GRID.Ny,[]);

u_tildeprime = 100*reshape(u_tildeprime(:,1,:),GRID.Nx,GRID.Ny,[]);
v_tildeprime = 100*reshape(v_tildeprime(:,1,:),GRID.Nx,GRID.Ny,[]);

ubar = 100*reshape(ubar,GRID.Nx,GRID.Ny,[]);
vbar = 100*reshape(vbar,GRID.Nx,GRID.Ny,[]);

Tdepth = reshape(T(:,indx_pycnocline(1),:),GRID.Nx,GRID.Ny,[]);
zdepth = z(indx_pycnocline(1));
T = reshape(T(:,1,:),GRID.Nx,GRID.Ny,[]);



xv = reshape(xv,GRID.Nx,GRID.Ny)/1000;
yv = reshape(yv,GRID.Nx,GRID.Ny)/1000;
x_sponge = x_sponge/1000;

starttime = getvalue([datadir,'/suntans.dat'],'starttime');
mtime = t/24/3600+datenum(sprintf('%.6f',starttime),'yyyymmdd.HHMMSS');

% find midpoint in y

save('./SUNTANS_results_3Dmovie.mat','-v7.3');