function [] = SUNTANS_plot_quad()
close all

load SUNTANS_quad_results.mat

%% top down plot
load jet_white.mat

figure(3124)

% subplot(2,1,1)
% z=Quad.eta(:,:,end);
% c1 = [min(min(z)) max(max(z))];
% pcolorjw(Quad.x,Quad.y,z)
% hold on
% contour(Quad.x,Quad.y,Quad.depth,10,'color',0.5*[1 1 1]);
% caxis(c1);
% cb = colorbar;
% ylabel(cb,'\eta')
% axis equal

L = [0.9 0.8 0.7 0.5 0.15];

for i=1:length(L)
subplot(5,1,i)
% qsx = 25;
% qsy = 10;

u = Quad.u(:,:,round(L(i)*end),end);
v = Quad.v(:,:,round(L(i)*end),end);

qsx = max([1 round(size(u,1)/30)]);
qsy = max([1 round(size(u,2)/10)]);


z=sqrt(u.^2+v.^2);
% z=-Quad.depth(:,:);
c1 = [min(min(z)) max(max(z))];
if c1(1)==c1(2)
    c1(2) =c1(2)+1E-10;
end
pcolorjw(Quad.x,Quad.y,z)
hold on
contour(Quad.x,Quad.y,Quad.depth,10,'color',0.5*[1 1 1]);
hold on
x = Quad.x(1:qsx:end,1:qsy:end);
y = Quad.y(1:qsx:end,1:qsy:end);
u = Quad.u(1:qsx:end,1:qsy:end,round(L(i)*end),end);
v = Quad.v(1:qsx:end,1:qsy:end,round(L(i)*end),end);
quiver(x,y,u,v,0.5,'color','k')
caxis(c1);
colormap(jet)
cb = colorbar;
axis equal

ylabel(cb,'u')
if i<length(L)
set(gca,'xtick',[])
end
ylabel('z (m)')

end
print -djpeg -r300 figure_quad_ubar

%% slice
load jet_white 

figure(3542)

L = [0.1 0.3 0.5 0.7 0.9];
c1 = reshape(Quad.u(:,:,:,end),[],1);
c1 = [quantile(c1,0.05) quantile(c1,0.95)];

for i=1:length(L)

indxy = round(size(Quad.y,2)*L(i));
if indxy<1
    indxy=1;
end

subplot(5,1,i)
x = Quad.x(:,indxy);
y = Quad.z;
x = x*ones(1,length(y));
y = ones(size(x,1),1)*y';
u = squeeze(Quad.u(:,indxy,:,end));
v = squeeze(Quad.w(:,indxy,:,end));
c1 = reshape(sqrt(u.^2+v.^2),[],1);
c1 = [quantile(c1,0.05) quantile(c1,0.95)];

pcolorjw(x,y,sqrt(u.^2+v.^2))
hold on
qsx = max([1 round(size(x,1)/30)]);
qsy = max([1 round(size(x,2)/10)]);


quiver(x(1:qsx:end,1:qsy:end),y(1:qsx:end,1:qsy:end),...
    u(1:qsx:end,1:qsy:end),v(1:qsx:end,1:qsy:end),0.2,'color','k')

caxis(c1)
colormap(jet)
cb = colorbar;
hold on

plot(Quad.x(:,indxy),-Quad.depth(:,indxy),'-k')


ylabel(cb,'|u|')
if i<length(L)
set(gca,'xtick',[])
end
ylabel('z (m)')

end
xlabel('x (m)')

print -djpeg -r300 figure_quad_slice

%% slice pressure
load jet_white 

figure(333542)

L = [0.1 0.3 0.5 0.7 0.9];
c1 = reshape(Quad.u(:,:,:,end),[],1);
c1 = [quantile(c1,0.05) quantile(c1,0.95)];

for i=1:length(L)

indxy = round(size(Quad.y,2)*L(i));
if indxy<1
    indxy=1;
end

subplot(5,1,i)
x = Quad.x(:,indxy);
y = Quad.z;
x = x*ones(1,length(y));
y = ones(size(x,1),1)*y';
u = squeeze(Quad.u(:,indxy,:,end));
v = squeeze(Quad.w(:,indxy,:,end));
z = squeeze(Quad.p_dyn(:,indxy,:,end));
c1 = reshape(z,[],1);
c1 = [quantile(c1,0.05) quantile(c1,0.95)];
if c1(1)>=c1(2)
    c1(2)=1.0001*abs(c1(1));
end
if isempty(c1)
    c1=[1 1];
end
pcolorjw(x,y,z)
hold on
qsx = max([1 round(size(x,1)/30)]);
qsy = max([1 round(size(x,2)/10)]);


quiver(x(1:qsx:end,1:qsy:end),y(1:qsx:end,1:qsy:end),...
    u(1:qsx:end,1:qsy:end),v(1:qsx:end,1:qsy:end),0.2,'color','k')

hold on
plot(Quad.x(:,indxy),-Quad.depth(:,indxy),'-k')

caxis(c1)
colormap(french)
cb = colorbar;

ylabel(cb,'p_{dyn}')
if i<length(L)
set(gca,'xtick',[])
end
ylabel('z (m)')

end
xlabel('x (m)')

print -djpeg -r300 figure_quad_slice_pressure

%% slice vorticity
load jet_white 

figure(333)

L = [0.1 0.3 0.5 0.7 0.9];
c1 = reshape(Quad.u(:,:,:,end),[],1);
c1 = [quantile(c1,0.05) quantile(c1,0.95)];

for i=1:length(L)

indxy = round(size(Quad.y,2)*L(i));
if indxy<1
    indxy=1;
end

subplot(5,1,i)
x = Quad.x(:,indxy);
y = Quad.z;
x = x*ones(1,length(y));
y = ones(size(x,1),1)*y';
u = squeeze(Quad.u(:,indxy,:,end));
v = squeeze(Quad.w(:,indxy,:,end));
z = squeeze(Quad.vorty_bar(:,indxy,:,end));
c1 = reshape(z,[],1);
c1 = [quantile(c1,0.05) quantile(c1,0.95)];

pcolorjw(x,y,z)
hold on
qsx = max([1 round(size(x,1)/30)]);
qsy = max([1 round(size(x,2)/10)]);


quiver(x(1:qsx:end,1:qsy:end),y(1:qsx:end,1:qsy:end),...
    u(1:qsx:end,1:qsy:end),v(1:qsx:end,1:qsy:end),0.2,'color','k')

hold on
plot(Quad.x(:,indxy),-Quad.depth(:,indxy),'-k')

caxis(c1)
colormap(french)
cb = colorbar;

ylabel(cb,'\omega_y')
if i<length(L)
set(gca,'xtick',[])
end
ylabel('z (m)')

end
xlabel('x (m)')

print -djpeg -r300 figure_quad_slice_vorticity

%% slice vorticity
load jet_white 

figure(33673)

L = [0.1 0.3 0.5 0.7 0.9];
c1 = reshape(Quad.u(:,:,:,end),[],1);
c1 = [quantile(c1,0.05) quantile(c1,0.95)];

for i=1:length(L)

indxy = round(size(Quad.y,2)*L(i));
if indxy<1
    indxy=1;
end

subplot(5,1,i)
x = Quad.x(:,indxy);
y = Quad.z;
x = x*ones(1,length(y));
y = ones(size(x,1),1)*y';
u = squeeze(Quad.u(:,indxy,:,end));
v = squeeze(Quad.w(:,indxy,:,end));
z = squeeze(Quad.nuTbar(:,indxy,:));
c1 = reshape(z,[],1);
c1 = [quantile(c1,0.05) quantile(c1,0.95)];

if c1(1)==c1(2)
    c1(2)=c1(1)*1.00001;
end
pcolorjw(x,y,z)
hold on
qsx = max([1 round(size(x,1)/30)]);
qsy = max([1 round(size(x,2)/10)]);


quiver(x(1:qsx:end,1:qsy:end),y(1:qsx:end,1:qsy:end),...
    u(1:qsx:end,1:qsy:end),v(1:qsx:end,1:qsy:end),0.2,'color','k')

hold on
plot(Quad.x(:,indxy),-Quad.depth(:,indxy),'-k')

caxis(c1)
colormap(french)
cb = colorbar;

ylabel(cb,'\nu_T')
if i<length(L)
set(gca,'xtick',[])
end
ylabel('z (m)')

end
xlabel('x (m)')

print -djpeg -r300 figure_quad_slice_nuT


%% drag parameters
figure(243)

subplot(4,1,1)
xplot = reshape(Quad.x,[],1);
yplot = reshape(Quad.eta(:,:,end),[],1);

plot(xplot,yplot,'.','color',0.5*[1 1 1])
hold on
plot(Quad.xfit,Quad.yfit,'-k')
text(0.1,0.2,['<Cd>(u b/c) = ' num2str(Quad.Cd_avg)],...
    'units','normalized')
ylabel('\eta')

subplot(4,1,2)
xplot = reshape(Quad.x,[],1);
yplot = reshape(Quad.U(:,:,end),[],1);
plot(xplot,yplot,'.','color',0.5*[1 1 1])
xplot = squeeze(nanmean(Quad.x(:,:),2));
yplot = squeeze(nanmean(Quad.U(:,:,end),2));
hold on
plot(xplot,yplot,'-k')
ylabel('U')

subplot(4,1,3)
xplot = reshape(Quad.x,[],1);
yplot = reshape(Quad.V(:,:,end),[],1);
plot(xplot,yplot,'.','color',0.5*[1 1 1])
xplot = squeeze(nanmean(Quad.x(:,:),2));
yplot = squeeze(nanmean(Quad.V(:,:,end),2));
hold on
plot(xplot,yplot,'-k')
ylabel('V')

subplot(4,1,4)
xplot = reshape(Quad.x,[],1);
yplot = sqrt(Quad.U(:,:,end).^2+Quad.V(:,:,end).^2);
yplot = yplot.*(Quad.depth+Quad.eta(:,:,end));
yplot = reshape(yplot,[],1);
plot(xplot,yplot,'.','color',0.5*[1 1 1])
xplot = squeeze(nanmean(Quad.x(:,:),2));
yplot = sqrt(Quad.U(:,:,end).^2+Quad.V(:,:,end).^2);
yplot = yplot.*(Quad.depth+Quad.eta(:,:,end));
yplot = squeeze(nanmean(yplot,2));
hold on
plot(xplot,yplot,'-k')
ylabel('|q|')
xlabel('x')

print -djpeg -r300 figure_quad_Cd


%%
figure(5324)
if Quad.t(end)>24*3600
    tplot = Quad.t/24/3600;
    xlbl = 't (days)';
else
    tplot = Quad.t;
    xlbl = 't (s)';
end
DRAG.Cd_bulk_avg(1)=0;
DRAG.Cd_taux_avg(1)=0;
DRAG.Cd_px_avg(1)=0;
subplot(3,1,1)
plot(tplot,DRAG.Cd_bulk_avg,'-k',...
    tplot,DRAG.Cd_taux_avg,'--r',...
    tplot,DRAG.Cd_px_avg,'-.b',...
    tplot,DRAG.Cd_px_avg+DRAG.Cd_taux_avg,'-.g')
ylabel('<C_D> periodic')
text(0.8,0.8,['<C_D>_{avg}=' num2str(DRAG.Cd_bulk_avg_bar)],...
    'units','normalized')
lg=legend('bulk','\tau','p','\tau+p');
set(lg,'location','northwest','orientation','horizontal','box','off')

subplot(3,1,2)
plot(tplot,Quad.U_avg,'-k')
hold on
plot(tplot,Quad.U_avg+Quad.U_std,'color',0.5*[1 1 1])
hold on
plot(tplot,Quad.U_avg-Quad.U_std,'color',0.5*[1 1 1])
ylabel('U (m/s)')

subplot(3,1,3)
plot(tplot,Quad.V_avg,'-k')
hold on
plot(tplot,Quad.V_avg+Quad.V_std,'color',0.5*[1 1 1])
hold on
plot(tplot,Quad.V_avg-Quad.V_std,'color',0.5*[1 1 1])
ylabel('V (m/s)')

xlabel(xlbl)

print -djpeg -r300 figure_quad_Cd_time

%%
figure(23423)

subplot(3,1,1)

z=DRAG.Cd_px_bar;
c1 = [quantile(reshape(z,[],1),0.05) quantile(reshape(z,[],1),0.95)];
if c1(1)==c1(2)
    c1(2) =c1(2)+1E-10;
end
pcolorjw(Quad.x,Quad.y,z)
hold on
contour(Quad.x,Quad.y,Quad.depth,5,'color',0.5*[1 1 1]);

caxis(c1);
colormap(jet)
cb = colorbar;
axis equal

ylabel(cb,'Cdp')
set(gca,'xtick',[])
ylabel('z (m)')

subplot(3,1,2)
z=DRAG.Cd_taux_bar;
c1 = [quantile(reshape(z,[],1),0.05) quantile(reshape(z,[],1),0.95)];

if c1(1)==c1(2)
    c1(2) =c1(2)+1E-10;
end
pcolorjw(Quad.x,Quad.y,z)
hold on
contour(Quad.x,Quad.y,Quad.depth,5,'color',0.5*[1 1 1]);

caxis(c1);
colormap(jet)
cb = colorbar;
axis equal

ylabel(cb,'Cd\tau')
set(gca,'xtick',[])
ylabel('z (m)')

subplot(3,1,3)
z=DRAG.Cd_bulk_bar;
c1 = [quantile(reshape(z,[],1),0.05) quantile(reshape(z,[],1),0.95)];

if c1(1)==c1(2)
    c1(2) =c1(2)+1E-10;
end
pcolorjw(Quad.x,Quad.y,z)
hold on
contour(Quad.x,Quad.y,Quad.depth,5,'color',0.5*[1 1 1]);

caxis(c1);
colormap(jet)
cb = colorbar;
axis equal

ylabel(cb,'Cdb')
set(gca,'xtick',[])
ylabel('z (m)')

print -djpeg -r300 figure_quad_Cd_space

%%

figure(423)

z=squeeze(Quad.pbot(:,:,end)+9.81*1000*Quad.eta(:,:,end));
c1 = [quantile(reshape(z,[],1),0.05) quantile(reshape(z,[],1),0.95)];
if c1(1)==c1(2)
    c1(2) =c1(2)+1E-10;
end
pcolorjw(Quad.x,Quad.y,z)
hold on
contour(Quad.x,Quad.y,Quad.depth,5,'color',0.5*[1 1 1]);

caxis(c1);
colormap(french)
cb = colorbar;
axis equal

ylabel(cb,'P_{dyn,b}')
set(gca,'xtick',[])
ylabel('z (m)')

print -djpeg -r300 figure_bottom_pressure

%% momentum terms
load jet_white

figure(4523454)

lims = max([Quad.USx_bar Quad.NL1_bar Quad.NL2_bar Quad.PGx_bar Quad.Fs_bar Quad.BT_bar]);

lims = lims*[-1 1];



for i=1:8

   if i==1
        z = Quad.USx;
        txt = 'US';
        tp ='-';
    elseif i==2
        z = Quad.NL1;
        txt = 'NL1';
        tp ='--';
    elseif i==3
        z = Quad.NL2;
        txt = 'NL2';
        tp ='-.';
    elseif i==4
        z = Quad.PGx;
        txt = 'PG';
        tp ='-';
    elseif i==5
        z = Quad.Fs;
        txt = 'Fs';
        tp ='--';
    elseif i==6
        z = Quad.BT;
        txt = 'BT';
        tp ='-.';
    elseif i==7
        z = Quad.Fpx;
        txt = 'Fp';
        tp ='-.';
    elseif i==8
        z = Quad.Err;
        txt = 'Err';
        tp ='-.';
   end

subplot(4,2,i)
pcolorjw(Quad.x,Quad.y,z)
colorbar
colormap(jet_white)
caxis(lims)
text(0.01,1.02,txt,'units','normalized')

end

print -djpeg -r300 figure_quad_2D_momentum

%% plot profile of momemtum terms

figure(6345)

subplot(2,1,1)
yindx = round(size(Quad.x,2)/2);

for i=1:10

   if i==1
        z = Quad.USx;
        txt = 'US';
        tp ='-';
    elseif i==2
        z = Quad.NL1;
        txt = 'NL1';
        tp ='--';
    elseif i==3
        z = Quad.NL2;
        txt = 'NL2';
        tp ='-.';
    elseif i==4
        z = Quad.PGx;
        txt = 'PG';
        tp ='-';
    elseif i==5
        z = Quad.Fs;
        txt = 'Fs';
        tp ='--';
    elseif i==6
        z = Quad.BT;
        txt = 'BT';
        tp ='-.';
    elseif i==7
        z = Quad.Fpx;
        txt = 'Fp';
        tp ='-.';
    elseif i==8
        z = Quad.HD;
        txt = 'HD';
        tp ='--';
    elseif i==9
        z = Quad.ReS;
        txt = 'ReS';
        tp ='-';    
    elseif i==10
        z = Quad.Err;
        txt = 'Err';
        tp ='-.';
   end
   
   subplot(3,1,1)
   plot(Quad.x(:,yindx),z(:,yindx,end),tp)
   hold on
   
   subplot(3,1,2)
   plot(Quad.x(:,yindx),nanmean(z(:,:,end),2),tp)
   hold on
   

end
subplot(3,1,1)
ylabel('m/s^2')
text(0.01,0.9,['y = ' num2str(Quad.y(1,yindx))],'units','normalized')

subplot(3,1,2)
ylabel('m/s^2')
text(0.01,0.9,['spatial avg'],'units','normalized')
lg=legend('US','NL1','NL2','PG','Fs','BT','Fp','HD','Res','Err');
set(lg,'location','best','orientation','horizontal','box','off')


subplot(3,1,3)
scatter(reshape(Quad.x,[],1),reshape(-Quad.depth,[],1),5,reshape(Quad.y,[],1))
cb=colorbar;
ylabel('y')
set(cb,'location','east')
hold on
plot(Quad.x(:,yindx),-Quad.depth(:,yindx))
ylabel('z (m)')
xlabel('x (m)')

print -djpeg -r300 figure_quad_2D_slice_momenum

%% profile plot
figure(576)
xlims = [min(DRAG.PROF.ubar-2*DRAG.PROF.ustd) max(DRAG.PROF.ubar+2*DRAG.PROF.ustd)];

semilogy(DRAG.PROF.ubar,DRAG.PROF.z,'-k')
hold on
semilogy(DRAG.PROF.ubar-2*DRAG.PROF.ustd,DRAG.PROF.z,'.','color',0.5*[1 1 1])

hold on
semilogy(xlims,DRAG.PROF.ztop*[1 1],'--','color',0.5*[1 1 1])
hold on
semilogy(xlims,DRAG.PROF.zmid*[1 1],'-.','color',0.5*[1 1 1])

hold on
semilogy(DRAG.PROF.ubar+2*DRAG.PROF.ustd,DRAG.PROF.z,'.','color',0.5*[1 1 1])

ylim([DRAG.PROF.z(end-1) DRAG.PROF.z(1)])

ylabel('z''(m)')
xlabel('u (m/s)')
lg=legend('<u>','+/-2\sigma_u','ztop','zmid');
set(lg,'location','northwest','box','on')

print -djpeg -r300 figure_quad_profile

%% spectral space of bottom
clear lgd 

figure(451635)

subplot(1,2,1)    
x = BATH.Sx;
y = BATH.Szx_bar;
if size(BATH.x,2)>3
    xx = BATH.Sy;
    yy = BATH.Szy_bar;
else
    xx = x;
    yy = nan+xx;
end
loglog(x,y,xx,yy)

text(0.1,0.9,['h_{rms}= ' num2str(BATH.hrms) ' m'],'units','normalized');


subplot(1,2,2)
x = BATH.Sx;
y = BATH.Sdzdx_bar;
if size(BATH.x,2)>3
    xx = BATH.Sy;
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
xlabel('k (m^{-1})')
% ylim([1e-7 1e-2])
text(0.02, 1.05,'(a) depth','units','normalized')

subplot(1,2,2)
ylabel('S_{dz_b/dx}')
xlabel('k (m^{-1})')
% ylim([1e-4 1e-2])
text(0.02, 1.05,'(b) slope','units','normalized')

legend('x','y','location','se')

% get(gcf,'position')
set(gcf,'position',[680   712   455   266])


print -djpeg -r300 figure_spectra_bottom

