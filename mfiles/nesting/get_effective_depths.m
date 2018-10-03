function get_effective_depths(h1,dmin,dmax)

% addpath '/Users/fringer/data/Research/normalmodes';

% km = 1000;
% hr = 3600;
h1=100;
dmin=1;
dmax=3000;
Nk = 500;
Ns = 20;
depths = linspace(dmin,dmax,Ns)';
c_s = zeros(Ns,1);
he_s = zeros(Ns,1);
for m=20%:Ns
    fprintf('On %d of %d\n',m,Ns);
    d0 = depths(m);
    dz = d0/Nk;
    zf = [0:dz:d0];
    
    g = 9.81;
    %    h1 = 75;
    h2 = d0-h1;
    delta_rho = 60;
    drho = 0.01;
    alpha_s = 0.99;
    
    rhof = -drho/2*tanh(2*atanh(alpha_s)/delta_rho*(zf-h2));
    drho_dz = diff(rhof)./diff(zf);
    z = 0.5*(zf(1:end-1)+zf(2:end));
    N = 1e-10+sqrt(-g*drho_dz);
    
    num_modes = 1;
    omega = 0;
    [phi,c1]=higher_modes(z(:),N(:),omega,num_modes,'rigid lid');
    c_s(m) = c1;
    phi = phi(:)/max(phi(:));
    z = z(:);
    
    dphi_dz = zeros(Nk,1);
    dphi_dz(2:Nk-1) = (phi(3:Nk)-phi(1:Nk-2))./(z(3:Nk)-z(1:Nk-2));
    dphi_dz(1) = 2*dphi_dz(2)-dphi_dz(3);
    dphi_dz(Nk) = 2*dphi_dz(Nk-1)-dphi_dz(Nk-2);
    
    he_s(m) = sqrt(3*trapz(z,phi.^2)/trapz(z,(dphi_dz).^2));
end

save heffective depths c_s he_s

