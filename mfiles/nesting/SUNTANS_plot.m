function [] = SUNTANS_plot()
disp('plotting results')
close all

load SUNTANS_results
GRID = load('./SUNTANS_grid.mat');

%%
set(0,'defaulttextinterpreter','latex');

% set(groot, 'defaultAxesTickLabelInterpreter','latex'); 

%% for 1d model, add some small random numbers to yv for plotting
indx = unique(yv);
if length(indx)==1
    yv = yv + 1e-5*randn(length(yv),1);
    ye = ye + 1e-5*randn(length(ye),1);
end
clear indx
%% try a profile plot
figure(55234)
clear uplot xplot yplot etaplot zbottom pplot wplot splot rhopplot rhobplot Uplot Wplot
load('./SUNTANS_results.mat','xv','x_sponge','x_wavemaker','Lx','sponge_distance');

load jet_white

i=round(Nout_keep);
j=round(GRID.Ny/2);
qsx = 20;
qsz = 4;

xplot = reshape(xv,GRID.Nx,GRID.Ny);
xplot = xplot(:,j);
[Zplot,Xplot]=meshgrid(z,xplot);

x_sponge = [x_sponge(1)+min(xplot), -x_sponge(1)+max(xplot)];

zbottom = -reshape(depth,GRID.Nx,GRID.Ny);
zbottom = zbottom(:,j);

etaplot = reshape(eta(:,1,i),GRID.Nx,GRID.Ny);
etaplot = etaplot(:,j);

plot1 = reshape(rho_b(:,:,i)+rho0,GRID.Nx,GRID.Ny,Nkmax);
plot1 = squeeze(plot1(:,j,:));

plot2 = reshape(rho_prime(:,:,i),GRID.Nx,GRID.Ny,Nkmax);
plot2 = squeeze(plot2(:,j,:));

uplot = reshape(u_tildeprime(:,:,i),GRID.Nx,GRID.Ny,Nkmax);
uplot = squeeze(uplot(:,j,:));

wplot = 0.5*(w(:,1:end-1,i)+w(:,2:end,i));
wplot = reshape(wplot,GRID.Nx,GRID.Ny,Nkmax);
wplot = squeeze(wplot(:,j,:));

Uplot = reshape(U(:,i),GRID.Nx,GRID.Ny);
Uplot = squeeze(Uplot(:,j))*ones(1,Nkmax);
Wplot = 0*Uplot;

Uplot(isnan(uplot))=nan;  

if max(xplot)-min(xplot)>2e3
    Xplot = Xplot/1000;
    xplot = xplot/1000;
    x_sponge = x_sponge/1000;
    x_wmkr = x_wavemaker/1000;
    xtxt = '(km)';
else
    xtxt = '(m)';
end
    
subplot(4,1,1)
cc = reshape(plot1,[],1);
cc = nanstd(cc);
pcolorjw(Xplot,Zplot,plot1)
colormap(jet_white)
cb=colorbar;
hold on
plot(xplot,etaplot,'-k')
hold on
plot(xplot,zbottom,'-k')
hold on
plot(x_sponge(1)*[1 1],[min(z) max(z)],'--k')
hold on
plot(x_sponge(2)*[1 1],[min(z) max(z)],'--k')
hold on
ylim([min(z)-0.05 0.05])
ylabel(cb,'\rho_b+\rho_0')
% xlabel('x ')
ylabel('z ')
title(['t= ' sprintf('%6.2f',t(i)/3600) ' hr'])
set(gca ,'Layer', 'Top')
set(gca,'xticklabel',[]);

subplot(4,1,2)
cc = reshape(plot2,[],1);
cc = nanstd(cc);
pcolorjw(Xplot,Zplot,plot2)
colormap(jet_white)
caxis(2*cc*[-1 1])
cb=colorbar;
hold on
contour(Xplot,Zplot,plot2,5,'color',0.5*[1 1 1])
hold on
plot(xplot,etaplot,'-k')
hold on
plot(xplot,zbottom,'-k')
hold on
plot(x_sponge(1)*[1 1],[min(z) max(z)],'--k')
hold on
plot(x_sponge(2)*[1 1],[min(z) max(z)],'--k')
hold on
ylim([min(z)-0.05 0.05])
ylabel(cb,'$\tilde{\rho}$','interpreter','latex')
ylabel('z ')
set(gca ,'Layer', 'Top')
set(gca,'xticklabel',[]);

subplot(4,1,3)
cc = reshape(Uplot,[],1);
cc = nanstd(cc);
pcolorjw(Xplot,Zplot,Uplot)
colormap(jet_white)
cb=colorbar;
hold on
plot(xplot,etaplot,'-k')
hold on
plot(xplot,zbottom,'-k')
hold on
plot(x_sponge(1)*[1 1],[min(z) max(z)],'--k')
hold on
plot(x_sponge(2)*[1 1],[min(z) max(z)],'--k')
hold on
ylim([min(z)-0.05 0.05])
ylabel(cb,'$\overline{U}+U''$','Interpreter','latex')
% xlabel('x/Lx ')
ylabel('z ')
set(gca ,'Layer', 'Top')
set(gca,'xticklabel',[]);

subplot(4,1,4)
cc = reshape(uplot,[],1);
cc = nanstd(cc);
pcolorjw(Xplot,Zplot,uplot)
colormap(jet_white)
caxis(2*cc*[-1 1])
cb=colorbar;
hold on
contour(Xplot,Zplot,uplot,5,'color',0.5*[1 1 1])
hold on
plot(xplot,etaplot,'-k')
hold on
plot(xplot,zbottom,'-k')
hold on
plot(x_sponge(1)*[1 1],[min(z) max(z)],'--k')
hold on
plot(x_sponge(2)*[1 1],[min(z) max(z)],'--k')
hold on
ylim([min(z)-0.05 0.05])
ylabel(cb,'$\tilde{u}''$','interpreter','latex')
xlabel(['x ' xtxt])
ylabel('z ')
set(gca ,'Layer', 'Top')

print -djpeg -r300 figure_slice

% zoom to pycnocline
for i=1:4
    subplot(4,1,i)
ylim([-D_pycnocline(3) 0])
end
print -djpeg -r300 figure_slice_zoom


close
%% slice of averages
% try a profile plot
if exist('AVG','var')
figure(54)
clear uplot uplot3 plot1 plot2 uplot_lp xplot yplot etaplot zbottom pplot wplot splot rhopplot rhobplot Uplot Wplot
load('./SUNTANS_results.mat','xv','x_sponge','x_wavemaker','Lx','sponge_distance');

load jet_white


i=round(AVG.Ntout);
j=round(GRID.Ny/2);
qsx = 20;
qsz = 4;

xplot = reshape(xv,GRID.Nx,GRID.Ny);
xplot = xplot(:,j);
[Zplot,Xplot]=meshgrid(z,xplot);

x_sponge = [x_sponge(1)+min(xplot), -x_sponge(1)+max(xplot)];

zbottom = -reshape(depth,GRID.Nx,GRID.Ny);
zbottom = zbottom(:,j);

etaplot = reshape(AVG.eta(:,1,i),GRID.Nx,GRID.Ny);
etaplot = etaplot(:,j);

plot1 = reshape(AVG.rho_b(:,:,i)+rho0,GRID.Nx,GRID.Ny,Nkmax);
plot1 = squeeze(plot1(:,j,:));

plot2 = reshape(AVG.rho_prime(:,:,i),GRID.Nx,GRID.Ny,Nkmax);
plot2 = squeeze(plot2(:,j,:));

uplot = reshape(AVG.u_prime(:,:,i),GRID.Nx,GRID.Ny,Nkmax);
uplot = squeeze(uplot(:,j,:));


wplot = 0.5*(AVG.w(:,1:end-1,i)+AVG.w(:,2:end,i));
wplot = reshape(wplot,GRID.Nx,GRID.Ny,Nkmax);
wplot = squeeze(wplot(:,j,:));

Uplot = reshape(AVG.U(:,i),GRID.Nx,GRID.Ny);
Uplot = squeeze(Uplot(:,j))*ones(1,Nkmax);
Wplot = 0*Uplot;

Uplot(isnan(uplot))=nan;  

if max(xplot)-min(xplot)>2e3
    Xplot = Xplot/1000;
    xplot = xplot/1000;
    x_sponge = x_sponge/1000;
    x_wmkr = x_wavemaker/1000;
    xtxt = '(km)';
else
    xtxt = '(m)';
end
    

subplot(3,1,1)
cc = reshape(plot2,[],1);
cc = nanstd(cc);
pcolorjw(Xplot,Zplot,plot2)
colormap(jet_white)
caxis(2*cc*[-1 1])
cb=colorbar;
hold on
contour(Xplot,Zplot,plot2,5,'color',0.5*[1 1 1])
hold on
plot(xplot,etaplot,'-k')
hold on
plot(xplot,zbottom,'-k')
hold on
plot(x_sponge(1)*[1 1],[min(z) max(z)],'--k')
hold on
plot(x_sponge(2)*[1 1],[min(z) max(z)],'--k')
hold on
ylim([min(z)-0.05 0.05])
ylabel(cb,'AVG $\tilde{\rho}$','interpreter','latex')
ylabel('z ')
set(gca ,'Layer', 'Top')
set(gca,'xticklabel',[]);

subplot(3,1,2)
cc = reshape(Uplot,[],1);
cc = nanstd(cc);
pcolorjw(Xplot,Zplot,Uplot)
colormap(jet_white)
cb=colorbar;
hold on
plot(xplot,etaplot,'-k')
hold on
plot(xplot,zbottom,'-k')
hold on
plot(x_sponge(1)*[1 1],[min(z) max(z)],'--k')
hold on
plot(x_sponge(2)*[1 1],[min(z) max(z)],'--k')
hold on
ylim([min(z)-0.05 0.05])
ylabel(cb,'AVG U')
% xlabel('x/Lx ')
ylabel('z ')
set(gca ,'Layer', 'Top')
set(gca,'xticklabel',[]);

subplot(3,1,3)
cc = reshape(uplot,[],1);
cc = nanstd(cc);
pcolorjw(Xplot,Zplot,uplot)
colormap(jet_white)
caxis(2*cc*[-1 1])
cb=colorbar;
hold on
contour(Xplot,Zplot,uplot,5,'color',0.5*[1 1 1])
hold on
plot(xplot,etaplot,'-k')
hold on
plot(xplot,zbottom,'-k')
hold on
plot(x_sponge(1)*[1 1],[min(z) max(z)],'--k')
hold on
plot(x_sponge(2)*[1 1],[min(z) max(z)],'--k')
hold on
ylim([min(z)-0.05 0.05])
ylabel(cb,'AVG $\tilde{u}$','interpreter','latex')
%xlabel(['x ' xtxt])
ylabel('z ')
set(gca ,'Layer', 'Top')


print -djpeg -r300 figure_slice_AVG

end
close
%% profile of velocities
figure(5633)
clear uplot xplot yplot etaplot zbottom pplot wplot splot rhopplot rhobplot Uplot Wplot plot1 plot2
load('./SUNTANS_results.mat','xv','x_sponge','x_wavemaker','Lx','sponge_distance');

load jet_white

i=round(Nout_keep);
j=round(GRID.Ny/2);
qsx = 20;
qsz = 4;

xplot = reshape(xv,GRID.Nx,GRID.Ny);
xplot = xplot(:,j);
[Zplot,Xplot]=meshgrid(z,xplot);

x_sponge = [x_sponge(1)+min(xplot), -x_sponge(1)+max(xplot)];

zbottom = -reshape(depth,GRID.Nx,GRID.Ny);
zbottom = zbottom(:,j);

etaplot = reshape(eta(:,1,i),GRID.Nx,GRID.Ny);
etaplot = etaplot(:,j);

plot1 = reshape(Ubar,GRID.Nx,GRID.Ny);
plot1 = squeeze(plot1(:,j))*ones(1,Nkmax);

plot2 = reshape(u_barprime,GRID.Nx,GRID.Ny,Nkmax);
plot2 = squeeze(plot2(:,j,:));

plot3 = reshape(Utilde(:,i),GRID.Nx,GRID.Ny);
plot3 = squeeze(plot3(:,j))*ones(1,Nkmax);

plot4 = reshape(u_tildeprime(:,:,i),GRID.Nx,GRID.Ny,Nkmax);
plot4 = squeeze(plot4(:,j,:));

uplot = reshape(u_tildeprime(:,:,i),GRID.Nx,GRID.Ny,Nkmax);
uplot = squeeze(uplot(:,j,:));

wplot = 0.5*(w(:,1:end-1,i)+w(:,2:end,i));
wplot = reshape(wplot,GRID.Nx,GRID.Ny,Nkmax);
wplot = squeeze(wplot(:,j,:));

Uplot = reshape(U(:,i),GRID.Nx,GRID.Ny);
Uplot = squeeze(Uplot(:,j))*ones(1,Nkmax);
Wplot = 0*Uplot;

Uplot(isnan(uplot))=nan;  

if max(xplot)-min(xplot)>2e3
    Xplot = Xplot/1000;
    xplot = xplot/1000;
    x_sponge = x_sponge/1000;
    x_wmkr = x_wavemaker/1000;
    xtxt = '(km)';
else
    xtxt = '(m)';
end
    
subplot(4,1,1)
plot1(isnan(plot2))=nan;
cc = reshape(plot1,[],1);
cc = nanstd(cc);
pcolorjw(Xplot,Zplot,plot1)
colormap(jet_white)
cb=colorbar;
hold on
plot(xplot,etaplot,'-k')
hold on
plot(xplot,zbottom,'-k')
hold on
plot(x_sponge(1)*[1 1],[min(z) max(z)],'--k')
hold on
plot(x_sponge(2)*[1 1],[min(z) max(z)],'--k')
hold on
ylim([min(z)-0.05 0.05])
ylabel(cb,'$\overline{U}$','interpreter','latex')
% xlabel('x ')
ylabel('z ')
title(['t= ' sprintf('%6.2f',t(i)/3600) ' hr'])
set(gca ,'Layer', 'Top')
set(gca,'xticklabel',[]);

subplot(4,1,2)
cc = reshape(plot2,[],1);
cc = nanstd(cc);
pcolorjw(Xplot,Zplot,plot2)
colormap(jet_white)
caxis(2*cc*[-1 1])
cb=colorbar;
hold on
contour(Xplot,Zplot,plot2,5,'color',0.5*[1 1 1])
hold on
plot(xplot,etaplot,'-k')
hold on
plot(xplot,zbottom,'-k')
hold on
plot(x_sponge(1)*[1 1],[min(z) max(z)],'--k')
hold on
plot(x_sponge(2)*[1 1],[min(z) max(z)],'--k')
hold on
ylim([min(z)-0.05 0.05])
ylabel(cb,'$\overline{u}''$','interpreter','latex')
ylabel('z ')
set(gca ,'Layer', 'Top')
set(gca,'xticklabel',[]);

subplot(4,1,3)
plot3(isnan(plot2))=nan;
cc = reshape(plot3,[],1);
cc = nanstd(cc);
pcolorjw(Xplot,Zplot,plot3)
colormap(jet_white)
cb=colorbar;
hold on
plot(xplot,etaplot,'-k')
hold on
plot(xplot,zbottom,'-k')
hold on
plot(x_sponge(1)*[1 1],[min(z) max(z)],'--k')
hold on
plot(x_sponge(2)*[1 1],[min(z) max(z)],'--k')
hold on
ylim([min(z)-0.05 0.05])
ylabel(cb,'$\tilde{U}$','interpreter','latex')
% xlabel('x/Lx ')
ylabel('z ')
set(gca ,'Layer', 'Top')
set(gca,'xticklabel',[]);

subplot(4,1,4)
cc = reshape(plot4,[],1);
cc = nanstd(cc);
pcolorjw(Xplot,Zplot,plot4)
colormap(jet_white)
caxis(2*cc*[-1 1])
cb=colorbar;
hold on
contour(Xplot,Zplot,plot4,5,'color',0.5*[1 1 1])
hold on
plot(xplot,etaplot,'-k')
hold on
plot(xplot,zbottom,'-k')
hold on
plot(x_sponge(1)*[1 1],[min(z) max(z)],'--k')
hold on
plot(x_sponge(2)*[1 1],[min(z) max(z)],'--k')
hold on
ylim([min(z)-0.05 0.05])
ylabel(cb,'$\tilde{u}''$','interpreter','latex')
xlabel(['x ' xtxt])
ylabel('z ')
set(gca ,'Layer', 'Top')

% end

print -djpeg -r300 figure_slice_velocities

close

%% density profile plot
figure(365234523)
Tplot = squeeze(nanmean(T(:,:,1),1));
Splot = squeeze(nanmean(S(:,:,1),1));
Rhoplot = squeeze(nanmean(rho(:,:,1),1));
Nplot = squeeze(nanmean(real(sqrt(N2(:,:,1))),1));

subplot(1,4,1)
plot(Tplot,z)
ylabel('z (m)')
xlabel('T (C)')

subplot(1,4,2)
plot(Splot,z)
xlabel('S (psu)')

subplot(1,4,3)
plot(Rhoplot,z)
xlabel('$\rho$ (kg/m$^3$)')

subplot(1,4,4)
plot(Nplot*3600/(2*pi),z)
xlabel('N (1/hr)')
print -djpeg -r300 figure_IC_profile
close
%% velocity profiles
figure(5334)
clear xplot indx
xplot(1) = min(xv);
xplot(2) = xplot(1)+sponge_distance/2;
xplot(3) = xplot(1)+1.5*sponge_distance;
xplot(4) = xplot(1)+0.33*Lx;
xplot(5) = xplot(1)+0.667*Lx;
xplot(6) = xplot(1)+Lx-1.5*sponge_distance;
xplot(7) = xplot(1)+Lx - 0.5*sponge_distance;
xplot(8) = xplot(1)+Lx;
ttl = {'L BC','L Spng','L wvmkr','0.33Lx','0.67Lx','R wvmkr','R Spng','R BC'};
N = length(xplot); % number of points
Nt = 2;
% xplot = linspace(min(xv),max(xv),N);
yplot = nanmean(yv);
for i=1:N % find closest point
   r = sqrt((xv-xplot(i)).^2 + (yv-yplot).^2);
   [~, indx(i)] = min(r);
end
indx2 = 1:Nt:length(t);
set(0,'DefaultAxesColorOrder',jet(length(indx2))) 

vmin = 0;
vmax = 0;
for i=1:N
    subplot(2,4,i)

        plot(squeeze(u(indx(i),:,indx2)),z)
        title(ttl{i});
%         title(['x=' num2str(xv(indx(i)))])
        hold on
        xlabel('u (m/s)')
        ylabel('z (m)')
        
        vmin = min([vmin nanmin(squeeze(u(indx(i),:,indx2)))]);
        vmax = max([vmax nanmax(squeeze(u(indx(i),:,indx2)))]);
end
% vmin = round(10*vmin)/10;
% vmax = round(10*vmax)/10;
if vmin==vmax
    vmax=vmin+1E-10;
end
for i=1:N
    subplot(2,4,i)
    xlim([vmin vmax])
end

print -djpeg -r300 figure_velocity_u_profile
close
%% velocity profiles
figure(5334)
clear xplot indx
xplot(1) = min(xv);
xplot(2) = xplot(1)+sponge_distance/2;
xplot(3) = xplot(1)+1.5*sponge_distance;
xplot(4) = xplot(1)+0.33*Lx;
xplot(5) = xplot(1)+0.667*Lx;
xplot(6) = xplot(1)+Lx-1.5*sponge_distance;
xplot(7) = xplot(1)+Lx - 0.5*sponge_distance;
xplot(8) = xplot(1)+Lx;
ttl = {'L BC','L Spng','L wvmkr','0.33Lx','0.67Lx','R wvmkr','R Spng','R BC'};
N = length(xplot); % number of points
Nt = 2;
% xplot = linspace(min(xv),max(xv),N);
yplot = nanmean(yv);
for i=1:N % find closest point
   r = sqrt((xv-xplot(i)).^2 + (yv-yplot).^2);
   [~, indx(i)] = min(r);
end
indx2 = 1:Nt:length(t);
set(0,'DefaultAxesColorOrder',jet(length(indx2))) 

vmin = 0;
vmax = 0;
for i=1:N
    subplot(2,4,i)

        plot(squeeze(v(indx(i),:,indx2)),z)
        title(ttl{i});
%         title(['x=' num2str(xv(indx(i)))])
        hold on
        xlabel('v (m/s)')
        ylabel('z (m)')
        
        vmin = min([vmin nanmin(squeeze(v(indx(i),:,indx2)))]);
        vmax = max([vmax nanmax(squeeze(v(indx(i),:,indx2)))]);
end
% vmin = round(10*vmin)/10;
% vmax = round(10*vmax)/10;
if vmin==vmax
    vmax=vmin+1E-10;
end
for i=1:N
    subplot(2,4,i)
    xlim([vmin vmax])
end

print -djpeg -r300 figure_velocity_v_profile
close

%% velocity profiles avg
if exist('AVG','var')
figure(539)
clear xplot indx
xplot(1) = min(xv);
xplot(2) = xplot(1)+sponge_distance/2;
xplot(3) = xplot(1)+1.5*sponge_distance;
xplot(4) = xplot(1)+0.33*Lx;
xplot(5) = xplot(1)+0.667*Lx;
xplot(6) = xplot(1)+Lx-1.5*sponge_distance;
xplot(7) = xplot(1)+Lx - 0.5*sponge_distance;
xplot(8) = xplot(1)+Lx;
ttl = {'L BC','L Spng','L wvmkr','0.33Lx','0.67Lx','R wvmkr','R Spng','R BC'};
N = length(xplot); % number of points
Nt = 2;
% xplot = linspace(min(xv),max(xv),N);
yplot = nanmean(yv);
for i=1:N % find closest point
   r = sqrt((xv-xplot(i)).^2 + (yv-yplot).^2);
   [~, indx(i)] = min(r);
end
indx2 = 1:length(AVG.t);
set(0,'DefaultAxesColorOrder',jet(length(indx2))) 

vmin = 0;
vmax = 0;
for i=1:N
    subplot(2,4,i)

        plot(squeeze(AVG.u(indx(i),:,indx2)),z)
        title(ttl{i});
%         title(['x=' num2str(xv(indx(i)))])
        hold on
        xlabel('AVG u (m/s)')
        ylabel('z (m)')
        
        vmin = min([vmin nanmin(squeeze(AVG.u(indx(i),:,indx2)))]);
        vmax = max([vmax nanmax(squeeze(AVG.u(indx(i),:,indx2)))]);
end
% vmin = round(10*vmin)/10;
% vmax = round(10*vmax)/10;
if vmin==vmax
    vmax=vmin+1E-10;
end
for i=1:N
    subplot(2,4,i)
    xlim([vmin vmax])
end

print -djpeg -r300 figure_velocity_profile_AVG
end
close
%% velocity profiles
figure(533488)
clear xplot indx
xplot(1) = min(xv);
xplot(2) = xplot(1)+sponge_distance/2;
xplot(3) = xplot(1)+1.5*sponge_distance;
xplot(4) = xplot(1)+0.33*Lx;
xplot(5) = xplot(1)+0.667*Lx;
xplot(6) = xplot(1)+Lx-1.5*sponge_distance;
xplot(7) = xplot(1)+Lx - 0.5*sponge_distance;
xplot(8) = xplot(1)+Lx;
ttl = {'L BC','L Spng','L wvmkr','0.33Lx','0.67Lx','R wvmkr','R Spng','R BC'};
N = length(xplot); % number of points
Nt = 2;
% xplot = linspace(min(xv),max(xv),N);
yplot = nanmean(yv);
for i=1:N % find closest point
   r = sqrt((xv-xplot(i)).^2 + (yv-yplot).^2);
   [~, indx(i)] = min(r);
end
indx2 = 1:size(u_tildeprime,3);
set(0,'DefaultAxesColorOrder',jet(length(indx2))) 

vmin = 0;
vmax = 0;
for i=1:N
    subplot(2,4,i)
    
%     plot(squeeze(0*ubar(indx(i),:)),z,'-k');
%     hold on

    qq = plot(squeeze(u_tildeprime(indx(i),:,indx2)),z);
    title(ttl{i});
%         title(['x=' num2str(xv(indx(i)))])
    hold on
%     qqq = plot(squeeze(ubar(indx(i),:)),z,'-k','linewidth',2);
    xlabel('$\tilde{u}''$ (m/s)')
    ylabel('$z$ (m)')

    vmin = min([vmin nanmin(squeeze(u_tildeprime(indx(i),:,indx2))) nanmin(squeeze(u_tildeprime(indx(i),:)))]);
    vmax = max([vmax nanmax(squeeze(u_tildeprime(indx(i),:,indx2))) nanmax(squeeze(u_tildeprime(indx(i),:)))]);
    grid on
end
% lg = legend([qq(1) qqq],'$u''$','$\overline{u}$');
% set(lg,'box','off','interpreter','latex')
% vmin = round(10*vmin)/10;
% vmax = round(10*vmax)/10;
if vmin==vmax
    vmax=vmin+1E-10;
end

for i=1:N
    subplot(2,4,i)
    xlim([vmin vmax])
end

print -djpeg -r300 figure_velocity_wave_profile
close
%% temperature profiles
figure(45236)
clear indx xplot yplot
xplot(1) = min(xv);
xplot(2) = xplot(1)+sponge_distance/2;
xplot(3) = xplot(1)+1.5*sponge_distance;
xplot(4) = xplot(1)+0.33*Lx;
xplot(5) = xplot(1)+0.667*Lx;
xplot(6) = xplot(1)+Lx-1.5*sponge_distance;
xplot(7) = xplot(1)+Lx - 0.5*sponge_distance;
xplot(8) = xplot(1)+Lx;
ttl = {'L BC','L Spng','L wvmkr','0.33Lx','0.67Lx','R wvmkr','R Spng','R BC'};
N = length(xplot); % number of points
Nt = 2;
% xplot = linspace(min(xv),max(xv),N);
yplot = nanmean(yv);
for i=1:N % find closest point
   r = sqrt((xv-xplot(i)).^2 + (yv-yplot).^2);
   [~, indx(i)] = min(r);
end
indx2 = 1:Nt:length(t);
set(0,'DefaultAxesColorOrder',jet(length(indx2))) 

for i=1:N
    subplot(2,4,i)

        plot(squeeze(T(indx(i),:,indx2)),z)
         title(ttl{i});
%         title(['x=' num2str(xv(indx(i)))])
        hold on
        xlabel('T (C)')
        ylabel('z (m)')
    
end

print -djpeg -r300 figure_temperature_profile
close
%% salt profiles
figure(54567)
clear indx xplot yplot

xplot(1) = min(xv);
xplot(2) = xplot(1)+sponge_distance/2;
xplot(3) = xplot(1)+1.5*sponge_distance;
xplot(4) = xplot(1)+0.33*Lx;
xplot(5) = xplot(1)+0.667*Lx;
xplot(6) = xplot(1)+Lx-1.5*sponge_distance;
xplot(7) = xplot(1)+Lx - 0.5*sponge_distance;
xplot(8) = xplot(1)+Lx;
ttl = {'L BC','L Spng','L wvmkr','0.33Lx','0.67Lx','R wvmkr','R Spng','R BC'};
N = length(xplot); % number of points
Nt = 2;
% xplot = linspace(min(xv),max(xv),N);
yplot = nanmean(yv);
for i=1:N % find closest point
   r = sqrt((xv-xplot(i)).^2 + (yv-yplot).^2);
   [~, indx(i)] = min(r);
end
indx2 = 1:Nt:length(t);
set(0,'DefaultAxesColorOrder',jet(length(indx2))) 


for i=1:N
    subplot(2,4,i)

        plot(squeeze(S(indx(i),:,indx2)),z)
        title(ttl{i});
%         title(['x=' num2str(xv(indx(i)))])
        hold on
        xlabel('S (psu)')
        ylabel('z (m)')
    
end

print -djpeg -r300 figure_salinty_profile
close
%% barotropic time
set(0,'DefaultAxesColorOrder','default') 
clear indx

figure(54324534)
xplot = [0 x_wavemaker(1)  Lx/2 x_wavemaker(2) Lx]+min(xv);
yplot = Ly/2+min(yv);
if t(end)>3600*24
    tday = t/3600/24;
    tday_low = t_lowpass/3600/24;
    xlbl = 't (day)';
else
    tday = t;
    tday_low=t_lowpass;
    xlbl = 't (s)';
end
for i=1:length(xplot)
   r = sqrt((xv-xplot(i)).^2 + (yv-yplot).^2);
   [~, indx(i)] = min(r);
end

subplot(6,1,1)
y = squeeze(eta(indx,:,:));
plot(tday,y)
lg = legend('L b/c', 'waveL','mid','waveR','R b/c');
set(lg,'box','off','location','best','orientation','horizontal')
hold on
ylabel('$\eta$(m)')
set(gca,'xticklabels',[]);

subplot(6,1,2)
y = squeeze(Ubar(indx,:));
plot(tday_low,y)
ylabel('$\overline{U}$ (m/s)','interpreter','latex')

set(gca,'xticklabels',[]);

subplot(6,1,3)
y = squeeze(Utilde(indx,:));
plot(tday,y)
ylabel('$\tilde{U}$(m/s)')
set(gca,'xticklabels',[]);

subplot(6,1,4)
y = squeeze(V(indx,:));
plot(tday,y)
ylabel('$V$(m/s)')
set(gca,'xticklabels',[]);

subplot(6,1,5)
y = squeeze(nanmean(T(indx,:,:),2));
plot(tday,y)
ylabel('T(C)')
set(gca,'xticklabels',[]);

subplot(6,1,6)
y = squeeze(nanmean(S(indx,:,:),2));
plot(tday,y)
ylabel('S(psu)')

% subplot(5,1,5)
% y = log10(squeeze(nanmean(abs(diff(eta(:,:,:))),1)));
% plot(tday,y)
% ylabel('$\log_{10}\overline{(\Delta\eta)}$')

xlabel(xlbl)
print -djpeg -r300 figure_time_barotropic
close
%% pycnocline time
set(0,'DefaultAxesColorOrder','default') 
clear indx

% pick a point midway between
indx_plot = round(mean(indx_pycnocline([1 2])));

figure(540056)
xplot = [x_wavemaker(1)  Lx/2 x_wavemaker(2)]+min(xv);
yplot = Ly/2+min(yv);
if t(end)>3600*24
    tday = t/3600/24;
    tday_avg = AVG.t/3600/24;
    xlbl = 't (day)';
else
    tday = t;
    tday_avg = AVG.t;
    xlbl = 't (s)';
end
for i=1:length(xplot)
   r = sqrt((xv-xplot(i)).^2 + (yv-yplot).^2);
   [~, indx(i)] = min(r);
end

subplot(6,1,1)
y = squeeze(u_tildeprime(indx,indx_plot,:));
plot(tday,y)
ylabel('$\tilde{u}''$ (m/s)','interpreter','latex')
lg = legend('waveL','mid','waveR');
set(lg,'box','off','location','best','orientation','horizontal')
title(['baroclinic in time, at z = ' num2str(z(indx_plot)) ' m'])
set(gca,'xticklabels',[]);

subplot(6,1,2)
y = squeeze(Utilde(indx,:));
plot(tday,y)
ylabel('$\tilde{U}$ (m/s)','interpreter','latex')
set(gca,'xticklabels',[]);

subplot(6,1,3)
y = squeeze(v(indx,indx_plot,:));
plot(tday,y)
ylabel('$v$ (m/s)','interpreter','latex')
set(gca,'xticklabels',[]);

subplot(6,1,4)
y = squeeze(w(indx,indx_plot,:));
plot(tday,y)
ylabel('w (m/s)')
set(gca,'xticklabels',[]);

subplot(6,1,5)
y = squeeze(T(indx,indx_plot,:));
plot(tday,y)
ylabel('T(C)')
set(gca,'xticklabels',[]);

subplot(6,1,6)
y = squeeze(S(indx,indx_plot,:));
plot(tday,y)
ylabel('S(psu)')

xlabel(xlbl)
print -djpeg -r300 figure_time_pycnocline
close
%% plot depth avg quantities
clear uplot eteplot
load('./SUNTANS_results.mat','xv','x_sponge','x_wavemaker','Lx','sponge_distance');

figure(8458)
xplot = linspace(min(xv),max(xv),GRID.Nx);
yplot = nanmean(yv)+0*xplot;
% [Zplot,Xplot]=meshgrid(z,xplot);
xlims = [min(xplot) max(xplot)];

data = squeeze(depth);
F = scatteredInterpolant(xv,yv,data,'linear');
zbottom = -F(xplot,yplot);


for i=Nout_keep
    data = squeeze(eta(:,1,i));
    F = scatteredInterpolant(xv,yv,data,'linear');
    etaplot = F(xplot,yplot);
    
    data = squeeze(z_pyn(:,i));
    F = scatteredInterpolant(xv,yv,data,'linear');
    zpyn_plot = F(xplot,yplot);
    
    data = squeeze(U(:,i));
    F = scatteredInterpolant(xv,yv,data,'linear');
    uplot = F(xplot,yplot);
    
    data = squeeze(V(:,i));
    F = scatteredInterpolant(xv,yv,data,'linear');
    vplot = F(xplot,yplot); 
    
    data = squeeze(u_tildeprime(:,1,i));
    F = scatteredInterpolant(xv,yv,data,'linear');
    utop = F(xplot,yplot);
    
    data = squeeze(u_tildeprime(:,round(0.25*end),i));
    F = scatteredInterpolant(xv,yv,data,'linear');
    u75 = F(xplot,yplot);
    
    data = squeeze(u_tildeprime(:,round(0.5*end),i));
    F = scatteredInterpolant(xv,yv,data,'linear');
    u50 = F(xplot,yplot);
    
    data = squeeze(u_tildeprime(:,round(0.75*end),i));
    F = scatteredInterpolant(xv,yv,data,'linear');
    u25 = F(xplot,yplot);
    
    data = squeeze(u_tildeprime(:,end,i));
    F = scatteredInterpolant(xv,yv,data,'linear');
    ubot = F(xplot,yplot);
    
    data = squeeze(u_barprime(:,1,end));
    F = scatteredInterpolant(xv,yv,data,'linear');
    utop2 = F(xplot,yplot);
    
    data = squeeze(u_barprime(:,round(0.25*end),end));
    F = scatteredInterpolant(xv,yv,data,'linear');
    u752 = F(xplot,yplot);
    
    data = squeeze(u_barprime(:,round(0.5*end),end));
    F = scatteredInterpolant(xv,yv,data,'linear');
    u502 = F(xplot,yplot);
    
    data = squeeze(u_barprime(:,round(0.75*end),end));
    F = scatteredInterpolant(xv,yv,data,'linear');
    u252 = F(xplot,yplot);
    
    data = squeeze(u_barprime(:,end,end));
    F = scatteredInterpolant(xv,yv,data,'linear');
    ubot2 = F(xplot,yplot);
end

if max(xplot)-min(xplot)>2e3
    xv = xv/1000;
    xplot = xplot/1000;
    x_sponge = x_sponge/1000;
    x_wmkr = x_wavemaker/1000;
    Lx = Lx/1000;
    sponge_distance = sponge_distance/1000;
    xtxt = '(km)';
    xlims=xlims/1000;
else
    xtxt = '(m)';
end

subplot(5,1,1)
plot(xv,squeeze(eta(:,:,end)),'.','color',0.5*[1 1 1])
hold on
yplot = etaplot;
plot(xplot,yplot,'-k')
hold on
plot(sponge_distance*[1 1],[min(yplot) max(yplot)],'--k')
hold on
plot((Lx-sponge_distance)*[1 1], [min(yplot) max(yplot)],'--k')
ylabel('$\eta$')
title(['T=' num2str(t(Nout_keep)/3600) ' hrs'])
xlim(xlims)

subplot(5,1,2)
yplot = zpyn_plot;
plot(xplot,yplot,'-k')
hold on
plot(sponge_distance*[1 1],[min(yplot) max(yplot)],'--k')
hold on
plot((Lx-sponge_distance)*[1 1], [min(yplot) max(yplot)],'--k')
ylabel('$ \zeta$')
% title(['T=' num2str(t(Nout_keep)/3600) ' hrs'])
xlim(xlims)

subplot(5,1,3)
yplot = [uplot;vplot];
plot(xplot,yplot)
hold on
plot(sponge_distance*[1 1],[min(min(yplot)) max(max(yplot))],'--k')
hold on
plot((Lx-sponge_distance)*[1 1], [min(min(yplot)) max(max(yplot))],'--k')
ylabel('U,V m/s')
lg=legend('U','V');
set(lg,'box','off','orientation','horizontal','location','best');
xlim(xlims)

subplot(5,1,4)
yplot=[utop; u75; u50; u25; ubot];
plot(xplot,yplot)
hold on
plot(sponge_distance*[1 1],[min(min(yplot)) max(max(yplot))],'--k')
hold on
plot((Lx-sponge_distance)*[1 1], [min(min(yplot)) max(max(yplot))],'--k')
ylabel('$\tilde{u}''$ m/s','interpreter','latex')
lg=legend('top','75','mid','25','bot');
set(lg,'box','off','orientation','horizontal','location','best');
xlim(xlims)

subplot(5,1,5)
yplot=[utop2; u752; u502; u252; ubot2];
plot(xplot,yplot)
hold on
plot(sponge_distance*[1 1],[min(min(yplot)) max(max(yplot))],'--k')
hold on
plot((Lx-sponge_distance)*[1 1], [min(min(yplot)) max(max(yplot))],'--k')
xlim(xlims)
ylabel('$\overline{u}''$ m/s','interpreter','latex')
% yplot = zbottom;
% plot(xplot,yplot)
% hold on
% plot(sponge_distance*[1 1],[min(yplot) max(yplot)],'--k')
% hold on
% plot((Lx-sponge_distance)*[1 1], [min(yplot) max(yplot)],'--k')
% ylabel('z (m)')
xlabel(['x ' xtxt])

print -djpeg -r300 figure_1D_results
close
%% how to plot 2D profile data

if ~isempty('PROF') && ~isempty(PROF)
  load jet_white

for i=1:length(PROF.dataNames) 
   loc=i; 
  xlims = [min(PROF.mtime(1,:)) max(PROF.mtime(1,:))];
  ylims = [ max([-PROF.d0(loc) -1000]) 0];
  
  
  figure(153324);
  subplot(5,1,1)
  plot(PROF.mtime(1,:),PROF.etaplot{i})
  xlim(xlims)
  datetick('x','keeplimits')
  set(gca,'xticklabels',[]);
  set(gca,'layer','top')
  ylabel('$\eta$ (m)')
  title(sprintf('Location %d, %s (D = %.0fm, x = %.0f, y = %.0f km)',...
      loc,PROF.dataNames{i},PROF.d0(i),PROF.dataX(loc)/1000,PROF.dataY(loc)/1000),...
      'interpreter','none');
  
  subplot(5,1,2)  
  pcolorjw(PROF.mtime,PROF.Z,...
      PROF.uplot{loc});
  ylim(ylims);
  xlim(xlims)
  datetick('x','keeplimits')
  set(gca,'xticklabels',[]);
  set(gca,'layer','top')
  ylabel('z (m)')
  cb=colorbar;
  ylabel(cb,'u (m/s)')
  colormap(jet_white)
  
  subplot(5,1,3)  
  pcolorjw(PROF.mtime,PROF.Z,...
      PROF.vplot{loc});
  ylim(ylims);
  xlim(xlims)
  datetick('x','keeplimits')
  set(gca,'xticklabels',[]);
  set(gca,'layer','top')
  ylabel('z (m)')
  cb=colorbar;
  ylabel(cb,'v (m/s)')
  
  subplot(5,1,4)  
  contourf(PROF.mtime,PROF.Z,...
      PROF.Tplot{loc},10);
  ylim(ylims);
  xlim(xlims)
  datetick('x','keeplimits')
  set(gca,'xticklabels',[]);
  set(gca,'layer','top')
  ylabel('z (m)')
  cb=colorbar;
  ylabel(cb,'T (C)')
   
  subplot(5,1,5)  
  contourf(PROF.mtime,PROF.Z,...
      PROF.Splot{loc},10);
  ylim(ylims);
  xlim(xlims)
  datetick('x','keeplimits')
  set(gca,'layer','top')
  ylabel('z (m)')
  cb=colorbar;
  ylabel(cb,'s (psu)')
  
  h2 = get(gca,'position');
  subplot(5,1,1)
  h1 = get(gca,'position');
  h1([3,4])=h2([3 4]);
  set(gca,'position',h1);
  
  txt = ['figure_profile_' PROF.dataNames{i}];
  print('-djpeg','-r300',txt);
 close
end



end

%% 2D barotropic  profiles
if exist('PROF') && ~isempty(PROF) && ~isempty(PROF.etaplot)
figure(1097856);
load jet_white


subplot(3,1,1)
for loc=1:length(PROF.u0)
  plot(PROF.mtime(1,:),...
      PROF.etaplot{loc},'-');
  hold on
  
  
  
end
datetickzoom('x')
set(gca,'xticklabels',[])
lg=legend(PROF.dataNames);
set(lg,'box','off','orientation','horizontal','location','best','interpreter','none')
ylabel('$\eta$ (m)')

subplot(3,1,2)
for loc=1:length(PROF.u0)
  plot(PROF.mtime(1,:),...
      nanmean(PROF.uplot{loc},1));
  hold on
  leg{loc} = sprintf('x = %.0f km',PROF.dataX(loc)/1000);
  shading flat;
end
datetickzoom('x')
set(gca,'xticklabels',[])
ylabel('U (m/s)')


subplot(3,1,3)
for loc=1:length(PROF.u0)
  plot(PROF.mtime(1,:),...
      nanmean(PROF.vplot{loc},1));
  hold on
  leg{loc} = sprintf('x = %.0f km',PROF.dataX(loc)/1000);
  shading flat;
end
datetickzoom('x')
ylabel('V (m/s)')

print -djpeg -r300 figure_profile_barotropic
close
end

%% plot energy flux
clear uplot Fx0_plot Fxp_plot etaplot zplot
load('./SUNTANS_results.mat','xv','x_sponge','x_wavemaker','Lx','sponge_distance');
set(0,'DefaultAxesColorOrder','default') 

i=round(Nout_keep);
j=round(GRID.Ny/2);
qsx = 20;
qsz = 4;

xplot = reshape(xv,GRID.Nx,GRID.Ny);
xplot = xplot(:,j);
[Zplot,Xplot]=meshgrid(z,xplot);

x_sponge = [x_sponge(1)+min(xplot), -x_sponge(1)+max(xplot)];

zbottom = -reshape(depth,GRID.Nx,GRID.Ny);
zbottom = zbottom(:,j);

etaplot = reshape(eta(:,1,i),GRID.Nx,GRID.Ny);
etaplot = etaplot(:,j);

Fx0bar_plot = reshape(Fx_0bar,GRID.Nx,GRID.Ny);
Fx0bar_plot = squeeze(Fx0bar_plot(:,j));

Fxpbar_plot = reshape(Fx_0bar,GRID.Nx,GRID.Ny);
Fxpbar_plot = squeeze(Fxpbar_plot(:,j));

Fxpbar_pos_plot = reshape(Fx_primebar_pos,GRID.Nx,GRID.Ny);
Fxpbar_pos_plot = squeeze(Fxpbar_pos_plot(:,j));

Fxpbar_neg_plot = reshape(Fx_primebar_neg,GRID.Nx,GRID.Ny);
Fxpbar_neg_plot = squeeze(Fxpbar_neg_plot(:,j));
    
% c = reshape(u,[],1);
% c = nanstd(c);

Fx0_plot = reshape(Fx_0,GRID.Nx,GRID.Ny,Nout_keep);
Fx0_plot = squeeze(Fx0_plot(:,j,:));

Fxp_plot = reshape(Fx_prime,GRID.Nx,GRID.Ny,Nout_keep);
Fxp_plot = squeeze(Fxp_plot(:,j,:));


% for i=1:Nout_keep
% %     data = sqrt(squeeze(Fx_0(:,i)).^2+squeeze(Fy_0(:,i)).^2);
%     data = squeeze(Fx_0(:,i));
%     F = scatteredInterpolant(xv,yv,data,'linear');
%     Fx0_plot(:,i) = F(xplot,yplot);
%     
% %     data = sqrt(squeeze(Fx_prime(:,i)).^2+squeeze(Fy_prime(:,i)).^2);
%     data = squeeze(Fx_prime(:,i));
%     F = scatteredInterpolant(xv,yv,data,'linear');
%     Fxp_plot(:,i) = F(xplot,yplot);
%     
%     data = squeeze(eta(:,1,i));
%     F = scatteredInterpolant(xv,yv,data,'linear');
%     etaplot = F(xplot,yplot);
% end
% for j=1:Nkmax % for each z level
%     data = squeeze(u_tilde(:,j,end));
%     F = scatteredInterpolant(xv,yv,data,'linear');
%     uplot(:,j) = F(xplot,yplot);
% end

if max(xplot)-min(xplot)>2e3
    Xplot = Xplot/1000;
    xplot = xplot/1000;
    x_sponge = x_sponge/1000;
    x_wmkr = x_wavemaker/1000;
    xtxt = '(km)';
else
    xtxt = '(m)';
end
xlims = [min(xplot) max(xplot)];

%
figure(23774)


subplot(4,1,1)
ax=plotyy(xplot,Fx0bar_plot,...
    xplot,Fxpbar_plot);
hold on
plot(xplot,0*xplot,'-k')
set(ax,'ylimmode','manual')
hold on
plot(x_sponge(1)*[1 1],1e8*[-1 1],'--k')
hold on
plot(x_sponge(2)*[1 1],1e8*[-1 1],'--k')
hold on
ylabel(ax(1),'$\bar{F_{0x}}$','interpreter','latex')
ylabel(ax(2),'$\bar{F_{x}^{''}}$','interpreter','latex')
xlim(ax(1),xlims)
xlim(ax(2),xlims)

subplot(4,1,2)
plot(xplot,Fxpbar_pos_plot,...
    xplot,Fxpbar_neg_plot,...
    xplot,Fxpbar_plot);
hold on
plot(xplot,0*xplot,'-k')
axis manual
hold on
plot(x_sponge(1)*[1 1],1e8*[-1 1],'--k')
hold on
plot(x_sponge(2)*[1 1],1e8*[-1 1],'--k')
hold on
ylabel('$\bar{F_{x}^{''}}$','interpreter','latex')
xlim(xlims)
lg=legend('+','-','sum');
set(lg,'orientation','horizontal','location','best','box','off')

subplot(4,1,3)
c = reshape(Fx0_plot,[],1);
c = nanstd(c);
load jet_white
colormap(jet_white)
pcolorjw(xplot,t/3600/12.42,Fx0_plot')
cb=colorbar;
set(cb,'location','east')
caxis(3*c*[-1 1])
set(ax,'ylimmode','manual')
hold on
plot(x_sponge(1)*[1 1],1e8*[-1 1],'--k')
hold on
plot(x_sponge(2)*[1 1],1e8*[-1 1],'--k')
hold on
hold on
plot(xlims,thetaramp/3600/12.42*[1 1],'--k')
hold on
plot(xlims,(t(end)-Tavg)/3600/12.42*[1 1],'--k')
xlim(xlims)
ylabel(cb,'Fx0')
ylabel('$t/T_{M2}$','interpreter','latex')

subplot(4,1,4)
c = reshape(Fxp_plot,[],1);
c = nanstd(c);
load jet_white
colormap(jet_white)
pcolorjw(xplot,t/3600/12.42,Fxp_plot')
cb=colorbar;
set(cb,'location','east')
caxis(1*c*[-1 1])
set(ax,'ylimmode','manual')
hold on
plot(x_sponge(1)*[1 1],1e8*[-1 1],'--k')
hold on
plot(x_sponge(2)*[1 1],1e8*[-1 1],'--k')
hold on
plot(xlims,thetaramp/3600/12.42*[1 1],'--k')
hold on
plot(xlims,(t(end)-Tavg)/3600/12.42*[1 1],'--k')
xlim(xlims)
ylabel(cb,'Fx''')
ylabel('$t/T_{M2}$','interpreter','latex')
xlabel(['$x$ ' xtxt],'interpreter','latex')

print -djpeg -r300 figure_energy_flux
close
%% wavemaker performance

% try a profile plot
% if wave_nesting==1
% figure(54)
% clear uplot uplot1 uplot2 uplot3 plot1 plot2 uplot_lp xplot yplot etaplot ...
%     zbottom pplot wplot splot rhopplot rhobplot Uplot Wplot u_tilde_LS
% load('./SUNTANS_results.mat','xv','x_sponge','x_wavemaker','Lx','sponge_distance');
% 
% load jet_white
% 
% 
% % theoretical wavemaker
% UHiw=0.50; Phiiw = 0;
% Tm2=12.42*3600;
% alpha1 = -0.7198; alpha2=-1.6810; delta=300; drho=5.5; h2=2500;
% fz = alpha1 * (1+alpha2-tanh((z+D-h2)/delta)).*dz;
% fzb = sum(fz)./sum(dz);
% fz = alpha1 * (1+alpha2-tanh((z+D-h2)/delta))-fzb;
% kiw = 2*pi/160e3;
% for i=1:size(AVG.uw_var_avg,1)
%        u_tilde_LS(i,:)=UHiw.*abs(fz).*abs(n1(i));
%    if n1(i)==0
%        u_tilde_LS(i,:)=UHiw.*abs(fz)+nan;
%    end
% end
% 
% % compute beta
% % u_tilde_rms = squeeze(2/sqrt(2)*sqrt(nanmean(u_tilde(:,:,61:81).^2,3)));
% u_tilde_rms = 2/sqrt(2)*sqrt(AVG.uw_var_avg(:,:,end));
% 
% beta = tanh(-u_tilde_rms./u_tilde_LS+1)+1;
% for i=1:size(beta,1)
%     beta_bar(i) = sum(beta(i,:).*dz')./sum(dz);
% end
% % u_LS_prime =  UHiw * fz * cos(kiw*xe-2*PI*prop->rtime/TM2+Phiiw);
% 
% %
% xplot = linspace(min(xv),max(xv),GRID.Nx);
% yplot = min(yv)+0.25*(max(yv)-min(yv))+0*xplot;
% [Zplot,Xplot]=meshgrid(z,xplot);
% qsx = 20;
% qsz = 4;
% data = squeeze(depth);
% F = scatteredInterpolant(xv,yv,data,'linear');
% zbottom = -F(xplot,yplot);
% 
% for i=round(AVG.Ntout)
%     data = squeeze(AVG.eta(:,1,i));
%     F = scatteredInterpolant(xv,yv,data,'linear');
%     etaplot = F(xplot,yplot);
%     
%     data = squeeze(beta_bar)';
%     F = scatteredInterpolant(xe,ye,data,'linear');
%     betaplot = F(xplot,yplot);
%     for j=1:Nkmax % for each z level
%         
%         data = rho(:,j,end);
%         indx44 = ~isnan(data);
%         data = data(indx44);
%         F = scatteredInterpolant(xv(indx44),yv(indx44),data,'linear');
%         uplot(:,j) = F(xplot,yplot);
%         
%         data = u_tilde_rms(:,j);
%         indx44 = ~isnan(data);
%         data = data(indx44);
%         F = scatteredInterpolant(xe(indx44),ye(indx44),data,'linear');
%         uplot1(:,j) = F(xplot,yplot);
%         
%         
%         data = u_tilde_LS(:,j);
%         indx44 = ~isnan(data);
%         data = data(indx44);
%         F = scatteredInterpolant(xe(indx44),ye(indx44),data,'linear');
%         uplot2(:,j) = F(xplot,yplot);
%                
%                
%         data = u_tilde_rms(:,j)./u_tilde_LS(:,j);
%         indx44 = ~isnan(data);
%         data = data(indx44);
%         F = scatteredInterpolant(xe(indx44),ye(indx44),data,'linear');
%         uplot3(:,j) = abs(F(xplot,yplot));
%         
%     end
% 
% 
% if max(xplot)-min(xplot)>2e3
%     Xplot = Xplot/1000;
%     xplot = xplot/1000;
%     x_sponge = x_sponge/1000;
%     x_wmkr = x_wavemaker/1000;
%     xtxt = '(km)';
% else
%     xtxt = '(m)';
% end
%     
% subplot(4,1,1)
% cc = reshape(uplot1,[],1);
% cc = nanstd(cc);
% pcolorjw(Xplot,Zplot,uplot1)
% colormap(jet_white)
% cb=colorbar;
% hold on
% plot(xplot,etaplot,'-k')
% hold on
% plot(xplot,zbottom,'-k')
% hold on
% plot(x_sponge(1)*[1 1],[min(z) max(z)],'--k')
% hold on
% plot(x_sponge(2)*[1 1],[min(z) max(z)],'--k')
% hold on
% plot(x_wmkr(1)*[1 1],[min(z) max(z)],'-.k')
% hold on
% plot(x_wmkr(2)*[1 1],[min(z) max(z)],'-.k')
% ylim([min(z)-0.05 0.05])
% ylabel(cb,'AVG $\tilde{u}''_{rms}$','interpreter','latex')
% % xlabel('x ')
% ylabel('z ')
% title(['Average NC, t= ' sprintf('%6.2f',AVG.t(i)/3600) ' hr'])
% set(gca ,'Layer', 'Top')
% set(gca,'xticklabel',[]);
% 
% subplot(4,1,2)
% cc = reshape(uplot2,[],1);
% cc = nanstd(cc);
% pcolorjw(Xplot,Zplot,uplot2)
% colormap(jet_white)
% % caxis(2*cc*[-1 1])
% cb=colorbar;
% hold on
% % contour(Xplot,Zplot,plot2,5,'color',0.5*[1 1 1])
% hold on
% plot(xplot,etaplot,'-k')
% hold on
% plot(xplot,zbottom,'-k')
% hold on
% plot(x_sponge(1)*[1 1],[min(z) max(z)],'--k')
% hold on
% plot(x_sponge(2)*[1 1],[min(z) max(z)],'--k')
% hold on
% plot(x_wmkr(1)*[1 1],[min(z) max(z)],'-.k')
% hold on
% plot(x_wmkr(2)*[1 1],[min(z) max(z)],'-.k')
% ylim([min(z)-0.05 0.05])
% ylabel(cb,'$\tilde{u}''_{rms,LS}$','interpreter','latex')
% ylabel('z ')
% set(gca ,'Layer', 'Top')
% set(gca,'xticklabel',[]);
% 
% subplot(4,1,3)
% cc = reshape(log10(uplot3),[],1);
% cc = nanstd(cc);
% pcolorjw(Xplot,Zplot,log10(uplot3))
% colormap(jet_white)
% cb=colorbar;
% caxis([-1 1])
% hold on
% plot(xplot,etaplot,'-k')
% hold on
% plot(xplot,zbottom,'-k')
% hold on
% plot(x_sponge(1)*[1 1],[min(z) max(z)],'--k')
% hold on
% plot(x_sponge(2)*[1 1],[min(z) max(z)],'--k')
% hold on
% plot(x_wmkr(1)*[1 1],[min(z) max(z)],'-.k')
% hold on
% plot(x_wmkr(2)*[1 1],[min(z) max(z)],'-.k')
% ylim([min(z)-0.05 0.05])
% ylabel(cb,'$\tilde{u}''_{rms}/\tilde{u}'',{rms,LS}$','interpreter','latex')
% % xlabel('x/Lx ')
% ylabel('z ')
% set(gca ,'Layer', 'Top')
% % set(gca,'xticklabel',[]);
% 
% 
% subplot(4,1,4)
% plot(xe/1000,beta_bar,'.k')
% hold on
% plot(xe/1000,beta(:,10),'.r')
% axis manual
% plot(x_sponge(1)*[1 1],[0 10],'--k')
% hold on
% plot(x_sponge(2)*[1 1],[0 10],'--k')
% hold on
% plot(x_wmkr(1)*[1 1],[0 10],'-.k')
% hold on
% plot(x_wmkr(2)*[1 1],[0 10],'-.k')
% ylabel('$\beta(x)$')
% lg=legend('$<\beta>$','$\beta(z=-100m)$');
% set(lg,'interpreter','latex')
% xlabel(['x ' xtxt])
% 
% 
% 
% drawnow
% pause(0.1)
% hold off
% end
% 
% print -djpeg -r300 figure_slice_wavemaker
% 
% end
% close

%% 2d plots
GRID = load('./SUNTANS_grid.mat');
if GRID.Ny>3 && GRID.Nx>3
xplot = reshape(xv,GRID.Nx,GRID.Ny)/1000;
yplot = reshape(yv,GRID.Nx,GRID.Ny)/1000;

uplot = reshape(u(:,1,end),GRID.Nx,GRID.Ny);
vplot = reshape(v(:,1,end),GRID.Nx,GRID.Ny);
Tplot = reshape(T(:,1,end),GRID.Nx,GRID.Ny);
Splot = reshape(S(:,1,end),GRID.Nx,GRID.Ny);
umag = sqrt(uplot.^2+vplot.^2);

zetaplot = reshape(z_pyn(:,end),GRID.Nx,GRID.Ny);
etaplot = reshape(eta(:,1,end),GRID.Nx,GRID.Ny);
qx = 10;

figure(23423)

subplot(2,2,1)
uplot = reshape(u_tildeprime(:,1,end),GRID.Nx,GRID.Ny);
vplot = reshape(v_tildeprime(:,1,end),GRID.Nx,GRID.Ny);
umag = sqrt(uplot.^2+vplot.^2);
pcolorjw(xplot,yplot,etaplot)
cb=colorbar;
colormap(dark_french)
ylabel(cb,'\eta(m)')
hold on
quiver(xplot(1:qx:end,1:qx:end),yplot(1:qx:end,1:qx:end),...
    uplot(1:qx:end,1:qx:end),vplot(1:qx:end,1:qx:end),0.5,'color','k');
ylabel('y (km)')
axis equal
text(0.01,1.1,'$\eta, u_{surf}$','units','normalized','interpreter','latex')

subplot(2,2,2)
uplot = reshape(U(:,end),GRID.Nx,GRID.Ny);
vplot = reshape(V(:,end),GRID.Nx,GRID.Ny);
umag = sqrt(uplot.^2+vplot.^2);
pcolorjw(xplot,yplot,umag)
cb=colorbar;
ylabel(cb,'|U| (m/s)')
hold on
quiver(xplot(1:qx:end,1:qx:end),yplot(1:qx:end,1:qx:end),...
    uplot(1:qx:end,1:qx:end),vplot(1:qx:end,1:qx:end),0.5,'color','k');
axis equal
text(0.01,1.1,'U','units','normalized')

subplot(2,2,3)
uplot = reshape(ubar(:,1),GRID.Nx,GRID.Ny);
vplot = reshape(vbar(:,1),GRID.Nx,GRID.Ny);
umag = sqrt(uplot.^2+vplot.^2);
pcolorjw(xplot,yplot,Tplot)
cb=colorbar;
ylabel(cb,'T (C)')
hold on
quiver(xplot(1:qx:end,1:qx:end),yplot(1:qx:end,1:qx:end),...
    uplot(1:qx:end,1:qx:end),vplot(1:qx:end,1:qx:end),0.5,'color','k');
ylabel('y (km)')
xlabel('x (km)')
axis equal
text(0.01,1.1,'Temp, $\overline{u}_{surf}$','units','normalized','interpreter','latex')

subplot(2,2,4)
uplot = reshape(ubar(:,1),GRID.Nx,GRID.Ny);
vplot = reshape(vbar(:,1),GRID.Nx,GRID.Ny);
umag = sqrt(uplot.^2+vplot.^2);
pcolorjw(xplot,yplot,Splot)
cb=colorbar;
ylabel(cb,'S (psu)')
hold on
quiver(xplot(1:qx:end,1:qx:end),yplot(1:qx:end,1:qx:end),...
    uplot(1:qx:end,1:qx:end),vplot(1:qx:end,1:qx:end),0.5,'color','k');
xlabel('x (km)')
axis equal
text(0.01,1.1,'Salt, $\overline{u}_{surf}$','units','normalized','interpreter','latex')

print -djpeg -r300 figure_2D_velocity_freesurface
close
end
%% plot the grid
% plot_grid(datadir)
close all

%%
disp('done plotting')
