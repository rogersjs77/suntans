
clear
close all

%
load SUNTANS_BC.mat

GRID = load('./SUNTANS_grid.mat');
%%

moviefile = 'movie_BC.mp4';
load jet_white
iskip=1;
qx=10;
cc = 0.3;
datadir='../data'; % data directory
sponge_distance = getvalue([datadir,'/suntans.dat'],'sponge_distance')/1000;
%%

tcount_all = size(BC.Fhat,3);
xv = BC.X/1000;
yv = BC.Y/1000;

xlimsBC = [min(min(xv)) max(max(xv))];
ylimsBC = [min(min(yv)) max(max(yv))];

xpt1 = GRID.BATH.corner_x/1000;
ypt1 = GRID.BATH.corner_y/1000;

xpt = GRID.BATH.sponge_x/1000;
ypt = GRID.BATH.sponge_y/1000;
% xlims = [min(BC.xv) max(BC.xv)]/1000;
% ylims = [min(BC.yv) max(BC.yv)]/1000;
% 
% xpt1 = [xlims(1),xlims(2),xlims(2),xlims(1),xlims(1)];
% ypt1 = [ylims(1),ylims(1),ylims(2),ylims(2),ylims(1)];
% 
% 
% xpt = [xlims(1)+sponge_distance,xlims(2)-sponge_distance,xlims(2)-sponge_distance,xlims(1)+sponge_distance,xlims(1)+sponge_distance];
% ypt = [ylims(1)+sponge_distance,ylims(1)+sponge_distance,ylims(2)-sponge_distance,ylims(2)-sponge_distance,ylims(1)+sponge_distance];

for i=1:iskip:tcount_all

Fx = BC.Fx(:,:,i);
Fy = BC.Fy(:,:,i);
F = sqrt(Fx.^2+Fy.^2);
Fhat = BC.Fhat(:,:,i);
eta = BC.eta(:,:,i);
% Fx = nanmean(BC.Fx(:,:,:),3);
% Fy = nanmean(BC.Fy(:,:,i),3);
% F = sqrt(Fx.^2+Fy.^2);
% Fhat = BC.Fhat(:,:,i);

subplot(2,1,1)

pcolorjw(xv,yv,eta)
hold on
quiver(xv(1:qx:end,1:qx:end),yv(1:qx:end,1:qx:end),...
    Fx(1:qx:end,1:qx:end),Fy(1:qx:end,1:qx:end),1,'k')
hold on
plot(xpt,ypt,'--','color',0.5*[1 1 1])
hold on
plot(xpt1,ypt1,'-','color',0.5*[1 1 1])
cb=colorbar;
colormap(dark_french)
ylabel(cb,'F''')
ylabel('y (km)')
axis equal
xlim(xlimsBC)
ylim(ylimsBC)
set(gca ,'Layer', 'Top')
title(['SUNTANS BC, F'', n = ' num2str(i), ', ' datestr(BC.mtime(i)) ])

F(Fhat==0)=nan;
subplot(2,1,2)

pcolorjw(xv,yv,F)
hold on
quiver(xv(1:qx:end,1:qx:end),yv(1:qx:end,1:qx:end),...
    Fx(1:qx:end,1:qx:end),Fy(1:qx:end,1:qx:end),1,'k')
hold on
plot(xpt,ypt,'--','color',0.5*[1 1 1])
hold on
plot(xpt1,ypt1,'-','color',0.5*[1 1 1])
cb=colorbar;
colormap(dark_french)
ylabel(cb,'F''Fhat')
ylabel('y (km)')
axis equal
xlim(xlimsBC)
ylim(ylimsBC)
set(gca ,'Layer', 'Top')
xlabel('x (km)')

print -djpeg -r150 movie.jpg
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

disp(['created ' num2str(i) ' of ' num2str(tcount_all)])

end
delete movie.jpg
close(vid);

%%


