function []=profile_points(datadir,name,lat,lon)
% compute xy points and save file
% Justin Rogers, 2018

% datadir = '../rundata';
% name = {'TC1','M2','M6','S7','B1'};
% lat = [21.0707 21.01 20.83 21.614 21.365];
% lon = [117.2202 118.16 119.05 117.282 118.594];

[N,E,Zone,lcm] = ell2utm(deg2rad(lat),deg2rad(lon));
if length(name)>length(N)
    GRID = load('./SUNTANS_grid.mat');
    offsetx = 0.01*GRID.L*cosd(GRID.theta);
    offsety = 0.01*GRID.L*sind(GRID.theta);
    indx = find(strcmp(name,'EBC'));
    if indx
        indx2 = 2;
        N(indx)=GRID.BATH.mid_y(indx2)-offsety;
        E(indx)=GRID.BATH.mid_x(indx2)-offsetx;       
    end
    indx = find(strcmp(name,'NBC'));
    if indx
        indx2 = 3;
        N(indx)=GRID.BATH.mid_y(indx2)-offsetx;
        E(indx)=GRID.BATH.mid_x(indx2)-offsety;      
    end
    indx = find(strcmp(name,'WBC'));
    if indx
        indx2 = 4;
        N(indx)=GRID.BATH.mid_y(indx2)+offsetx;
        E(indx)=GRID.BATH.mid_x(indx2)-offsety;    
    end
    indx = find(strcmp(name,'SBC'));
    if indx
        indx2 = 1;
        N(indx)=GRID.BATH.mid_y(indx2)+offsety;
        E(indx)=GRID.BATH.mid_x(indx2)+offsetx;      
    end
    if length(N)~=length(name)
        disp('ERROR: points not matching')
        return
    end
    
end



for i=1:length(N)
SUN_loc{i} = [ sprintf('%15.2f',E(i))...
                sprintf('%15.2f',N(i))]; 
end
saveascii(SUN_loc,[datadir '/dataxy.dat'])
saveascii(name,[datadir '/dataxy_names.dat'])

%% plot points on grid
load('./SUNTANS_grid.mat','xv','yv','Depth','BATH')
xv=xv/1000;
yv=yv/1000;
N=N/1000;
E=E/1000;

%% get sponge layer points
xpt = BATH.sponge_x/1000;
ypt = BATH.sponge_y/1000;
% x_sponge = getvalue([datadir,'/suntans.dat'],'sponge_distance')/1000;
% 
xlims = [min(min(xv)) max(max(xv))];
ylims = [min(min(yv)) max(max(yv))];
% 
% xpt = [xlims(1)+x_sponge(1),xlims(2)-x_sponge(1),xlims(2)-x_sponge(1),xlims(1)+x_sponge(1),xlims(1)+x_sponge(1)];
% ypt = [ylims(1)+x_sponge(1),ylims(1)+x_sponge(1),ylims(2)-x_sponge(1),ylims(2)-x_sponge(1),ylims(1)+x_sponge(1)];

%%
figure(234)
scatter(xv,yv,5,-Depth)
cb=colorbar;
ylabel(cb,'z (m)')
hold on
indx = Depth<30;
plot(xv(indx),yv(indx),'.k')
hold on
plot(xpt,ypt,'--','color',0*[1 1 1])
hold on
plot(E,N,'r*')
text(E+1.5,N,name,'interpreter','none')
xlabel('x (km UTM)')
ylabel('y (km UTM)')
axis equal
xlim(xlims)
ylim(ylims)

print -djpeg -r300 figure_points
close

end
