clear

load('./SUNTANS_results.mat',...
    'xv','yv','z','u_tildeprime','rho','z_pyn','indx_pycnocline','x_sponge','t','Nkmax');

GRID = load(['./SUNTANS_grid.mat']);

rho = reshape(rho,GRID.Nx,GRID.Ny,Nkmax,[]);
u_tildeprime = -100*reshape(u_tildeprime,GRID.Nx,GRID.Ny,Nkmax,[]);
xv = -(xv-max(xv));
yv = -(yv-max(yv));
xv = reshape(xv,GRID.Nx,GRID.Ny)/1000;
yv = reshape(yv,GRID.Nx,GRID.Ny)/1000;
z_pyn = reshape(z_pyn,GRID.Nx,GRID.Ny,[]);
x_sponge = x_sponge/1000;

% find midpoint in y
indxy = round(size(xv,2)/2);
x = xv(:,indxy);
[Z,X]=meshgrid(z,xv(:,indxy));

save('./SUNTANS_results_movie.mat','-v7.3');
