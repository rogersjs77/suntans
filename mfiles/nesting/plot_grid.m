function [] = plot_grid(datadir,txt)
close all

%% General Settings
% clear
% datadir='../rundata';
state = 'sal'; %'temp','sal','full'
EMPTY=999999; % code in SUNTANS data for uncomputed cells... do not change.
% datadir ='../../runs/3-9_h200_nonhs_witch/'; % directory of SUNTANS data
savename = 'SUNTANS_results'; % name of saved .mat file
savedir = './'; % directory to save data .mat file

makeplots = 0; % zero is off, 1 is on
rho0 = 1000;
grav = 9.81;

if nargin<2
    txt=0;
end

%% Bathymetry settings

depth=load([datadir,'/depth.dat']);
% % xv = depth(:,1);
% % yv = depth(:,2);
depth = depth(:,3);
% D = max(depth); %max depth


%% Loading grid and time data
c = load([datadir,'/cells.dat']);
if size(c,2)<9
    xv = c(:,1);
    yv = c(:,2);
else
    xv = c(:,2);
    yv = c(:,3);
end
points = load([datadir,'/points.dat']);
edges = load([datadir,'/edges.dat']);
% dz = load([datadir,'/vertspace.dat']);
% Nkmax = getvalue([datadir,'/suntans.dat'],'Nkmax');
% % z = getz(dz); %depth
% [Nc,~] = size(c); %[number of cells, number of columns in cells.dat


%% downsample
r=1;
% if length(xv) > 1E3
%    indx = round(1+(length(xv)-1)*rand(1E4,1));
%    r = round(length(xv)/length(indx));
%    xv = xv(indx);
%    yv = yv(indx);
%    
%    
% end


%%



figure(2342)

scatter(xv,yv,5,-depth')
cb = colorbar;
ylabel(cb,'z (m)')
% plot(xv,yv,'.k') %voronoi centers
hold on

if size(c,1)>1E3
    indx = find(edges(:,3)>0);
else
    indx = find(edges(:,3)>=0);
end

for ii=1:length(indx) % suntans uses C indexing starts at 0
   i=indx(ii);
   x = [points(edges(i,1)+1,1) points(edges(i,2)+1,1)];
   y = [points(edges(i,1)+1,2) points(edges(i,2)+1,2)];
   mk = edges(i,3);
   if mk==0
       cb = 'k';
   elseif mk==1
       cb = 'b';
   elseif mk==2
       cb = 'r';
   elseif mk==3
       cb = 'g';
   elseif mk==4
       cb = 'm';
   elseif mk==5
       cb = 'y';
   end
   plot(x,y,cb)
   hold on   
    
end

% plot(points(:,1),points(:,2),'ok');
axis equal
axis manual
hold on
p1= plot(-1E8*[1 1],-1E8*[1 1],'-k');hold on
p2=plot(-1E8*[1 1],-1E8*[1 1],'-b');hold on
p3=plot(-1E8*[1 1],-1E8*[1 1],'-r');hold on
p4=plot(-1E8*[1 1],-1E8*[1 1],'-g');hold on
p5=plot(-1E8*[1 1],-1E8*[1 1],'-m');hold on
ylabel('y (m)')
xlabel('x (m)')

qq=legend([p1 p2 p3 p4 p5],'0 comp','1 closed','2 open u','3 open \eta','4 no slip');
set(qq,'box','off','location','best')

if txt && size(c,1)<1E3
    for i =1:size(points,1)
    txt1{i} = num2str(i-1);
    end
    text(points(:,1)+10,points(:,2),txt1)
    
    clear txt1
    for i =1:size(c,1)
        txt1{i} = num2str(i-1);
    end
    if size(c,2)>8
        text(c(:,2)+10,c(:,3),txt1)
    else
        text(c(:,1)+10,c(:,2),txt1)
    end

    
end
% title(['model grid, every ' num2str(r) ' points shown'])

print -djpeg -r300 figure_grid_overview