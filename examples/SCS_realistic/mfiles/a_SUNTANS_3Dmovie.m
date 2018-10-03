
clear
close all
set(0,'defaulttextinterpreter','latex')
%
load SUNTANS_results_3Dmovie.mat
GRID = load('./SUNTANS_grid.mat');
%%

moviefile = 'movie_3D.mp4';
load jet_white
iskip=1;
qx=10;
cc = 0.3;

%%
tcount_all = length(t);

xlims = [min(min(xv)) max(max(xv))];
ylims = [min(min(yv)) max(max(yv))];

xpt = GRID.BATH.sponge_x/1000;
ypt = GRID.BATH.sponge_y/1000;
% xpt = [xlims(1)+x_sponge(1),xlims(2)-x_sponge(1),xlims(2)-x_sponge(1),xlims(1)+x_sponge(1),xlims(1)+x_sponge(1)];
% ypt = [ylims(1)+x_sponge(1),ylims(1)+x_sponge(1),ylims(2)-x_sponge(1),ylims(2)-x_sponge(1),ylims(1)+x_sponge(1)];

eta = eta - nanmean(reshape(eta,[],1));

etalims = nanmean(eta,3);
etalims = reshape(etalims,[],1);
etalims = [-1 1]*max(abs(etalims));
etac = linspace(etalims(1), etalims(2),8);

tlims = reshape(Tdepth,[],1);
tlims = [quantile(tlims,0.05) quantile(tlims,0.95)];
tc = linspace(tlims(1),tlims(2),8);

for i=1:iskip:tcount_all

ssh = eta(:,:,i);
ssh=ssh-nanmean(nanmean(ssh));
u = u_tildeprime(:,:,i);
v = v_tildeprime(:,:,i);

subplot(2,1,1)

pcolorjw(xv,yv,ssh)
hold on
contour(xv,yv,ssh,etac,'color',0.5*[1 1 1])
hold on
quiver(xv(1:qx:end,1:qx:end),yv(1:qx:end,1:qx:end),...
    u(1:qx:end,1:qx:end),v(1:qx:end,1:qx:end),1,'k')
hold on
plot(xpt,ypt,'--','color',0.5*[1 1 1])
cb=colorbar;
colormap(dark_french)
caxis(etalims)
% caxis(cc*[-1 1])
ylabel(cb,'$\eta''$ (m)','interpreter','latex')
ylabel('$y$ (km)')
axis equal
xlim(xlims)
ylim(ylims)
set(gca ,'Layer', 'Top')
title(['$\eta$, $\tilde{u}''(z = 0) , $' datestr(mtime(i),'dd-mmm-yyyy HH:MM')],'interpreter','latex')

ssT = Tdepth(:,:,i);
u = u_tildeprime_depth(:,:,i);
v = v_tildeprime_depth(:,:,i);

subplot(2,1,2)

pcolorjw(xv,yv,ssT)
hold on
contour(xv,yv,ssT,tc,'color',0.5*[1 1 1])
hold on
quiver(xv(1:qx:end,1:qx:end),yv(1:qx:end,1:qx:end),...
    u(1:qx:end,1:qx:end),v(1:qx:end,1:qx:end),1,'k')
hold on
plot(xpt,ypt,'--','color',0.5*[1 1 1])
cb=colorbar;
colormap(dark_french)
caxis(tlims)
ylabel(cb,'$T$ (C)','interpreter','latex')
ylabel('$y$ (km)')
axis equal
xlim(xlims)
ylim(ylims)
set(gca ,'Layer', 'Top')
title(['T, $\tilde{u}''(z = -140 $m) '],'interpreter','latex')
xlabel('$x$ (km)')

print -djpeg -r300 movie.jpg
close
% export_fig(['movie.jpg'],'-jpg','-transparent','-r140','-nocrop');

[im,map] = imread(['movie.jpg']);   
[imind,cm] = rgb2ind(im,256);
if i == 1
    vid = VideoWriter(moviefile,'MPEG-4');
    vid.FrameRate = 5;
    open(vid);
    writeVideo(vid,im);
%    imwrite(imind,cm,moviefile,'gif', 'Loopcount',inf);
else
%    imwrite(imind,cm,moviefile,'gif','WriteMode','append','delaytime',0.0);
   writeVideo(vid,im);
end      

disp(['created ' num2str(i) ' of ' num2str(length(t))])

end
delete movie.jpg
close(vid);

%%


