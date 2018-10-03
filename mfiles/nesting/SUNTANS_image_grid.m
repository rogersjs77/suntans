function [] = SUNTANS_image_grid(image_dir)


GRID = load('./SUNTANS_grid.mat');

[A,R]=geotiffread(image_dir);
info = geotiffinfo(image_dir);
[x_im_utm,y_im_utm]=pixcenters(info);

% format is Y,X, 
A = A(:,:,1:3);

% reduce data size to model area
indxi = x_im_utm > (min(reshape(GRID.BATH.x_UTM,[],1))) &...
    x_im_utm < (max(reshape(GRID.BATH.x_UTM,[],1)));
indxj = y_im_utm > (min(reshape(GRID.BATH.y_UTM,[],1))) &...
    y_im_utm < (max(reshape(GRID.BATH.y_UTM,[],1)));

%%
% figure(1)
% imagesc(x_im_utm,y_im_utm,A);
% set(gca,'ydir','normal')
% hold on
% plot(GRID.BATH.x_UTM(1:20:end),GRID.BATH.y_UTM(1:20:end),'.r')
% axis equal
% 
% figure(2)
% imagesc(x_im_utm(indxi),y_im_utm(indxj),A(indxj,indxi,:));
% set(gca,'ydir','normal')
% hold on
% plot(GRID.BATH.x_UTM(1:20:end),GRID.BATH.y_UTM(1:20:end),'.r')
% axis equal
%% reduce grid sizes for interpolation
A = A(indxj,indxi,:);
x_im_utm = x_im_utm(indxi);
y_im_utm = y_im_utm(indxj);
[X_im_utm,Y_im_utm] = meshgrid(x_im_utm,y_im_utm);
clear indxi indxj


%% format is Y,X
dx_im = x_im_utm(2)-x_im_utm(1);
L = GRID.L;
W = GRID.W;
theta = GRID.theta;

% rotate about midpoint
Ar = imrotate(A,-theta,'nearest');
% midpoint
indxi = round(size(Ar,1)/2);
indxj = round(size(Ar,2)/2);
% what points to keep
indxi = (indxi-round(W/dx_im/2)):1:(indxi+round(W/dx_im/2));
indxj = (indxj-round(L/dx_im/2)):1:(indxj+round(L/dx_im/2));

Im = Ar(indxi,indxj,:);
x = linspace(GRID.BATH.x(1,1),GRID.BATH.x(end,end),size(Ar,2));
y = linspace(GRID.BATH.y(end,end),GRID.BATH.y(1,1),size(Ar,1));

%%
figure(144)
subplot(2,1,1)
pcolorjw(GRID.BATH.x,GRID.BATH.y,GRID.BATH.Z)
cb = colorbar;
set(cb,'location','east');
xlabel(cb,'z (m)')
axis equal
colormap jet
ylabel('y (m)')

subplot(2,1,2)
imagesc(x,y,Im);
set(gca,'ydir','normal')
axis equal
ylabel('y (m)')
xlabel('x (m)')

print -djpeg -r300 figure_image_grid.jpg


%%
clear X_im Y_im A X_im_utm Y_im_utm A Ar cb GRID x_im_utm y_im_utm indxi indxj

%% save

save('./SUNTANS_image.mat','-v7.3','x','y','Im')
end