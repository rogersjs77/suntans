clear
% 04-Jun-2011 02:50:00
mtimesave = datenum(2011,6,4,2,50,0);

% load data in peices to minimize memory
load('./SUNTANS_results.mat',...
't','mtime','D','Nkmax','z','dz','N2','rho0','xv','yv','xe','ye',...
'eta','z_pyn','etabar','indx_pycnocline','Tavg','PROF',...
'Fx_0bar','Fx_primebar','Fx_0','Fx_prime','Fy_0bar','Fy_primebar',...
'x_sponge','x_wavemaker','depth','Nout_keep')


%% save only nearest time step

[~,indx] = min(abs(mtimesave-mtime));
mtime = mtime(indx);
tsave = t(indx);

N2_init = N2(:,:,1);
eta = eta(:,:,indx);
N2 = N2(:,:,indx);

load('./SUNTANS_results.mat','rho','rho_prime');
rho_init = rho(:,:,1);
rho = rho(:,:,indx);
rho_prime = rho_prime(:,:,indx);

load('./SUNTANS_results.mat','u','v');
u = u(:,:,indx);
v = v(:,:,indx);

load('./SUNTANS_results.mat',...
'Ubar','Utilde','u_tildeprime','u_barprime',...
'Vbar','Vtilde','v_tildeprime','v_barprime')
u_tildeprime = u_tildeprime(:,:,indx);
Utilde = Utilde(:,indx);
v_tildeprime = v_tildeprime(:,:,indx);
Vtilde = Vtilde(:,indx);
z_pyn = z_pyn(:,indx);


load('./SUNTANS_results.mat','T','S','bottom_cells','AVG')
S = S(:,:,indx);
for i=1:length(xv)
   Tbot(i,:) = squeeze(T(i,bottom_cells(i),:));
end
T = T(:,:,indx);
load('./SUNTANS_results.mat','AVG')
indx = AVG.t>(AVG.t(end)-Tavg);
etabar = squeeze(nanmean(AVG.eta(:,:,indx),3));

clear indx AVG i
save ./SUNTANS_results_reduce
