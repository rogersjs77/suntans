function [phiv] = interpNCOM_parallel(XN,YN,zN,phiN,xv,yv,zv)
% this function takes NCOM xN,yN,zN,phi data
% interpolates in xy using parfor to obtain phi (T,S,U,V)

% XN = NC.X;
% YN = NC.Y;
% zN = NC.z;
% phiN = NC.temp;
% zv=z_r;

% XN = reshape(XN,[],1);
% YN = reshape(YN,[],1);

if length(size(phiN))==3 %3d variable
    
    % interp in xy at each z level
    temp=[];
    phiv_z = [];
    parfor i=1:size(phiN,3)
        temp = squeeze(phiN(:,:,i));
        temp = reshape(temp,[],1);
        F = scatteredInterpolant(reshape(XN,[],1),reshape(YN,[],1),reshape(temp,[],1));
%         F = griddedInterpolant(XN,YN,temp,'linear');
        phiv_z(:,i) = F(xv,yv);   
    end
    
    % now interp in z, extrapolate unnknown values
    phiv = interp1(zN,phiv_z',zv,'linear','extrap')';
    
    phiv = naninterp(phiv,0);
    
elseif  length(size(phiN))==2 % 2d variable
    
        temp = squeeze(phiN);
        F = scatteredInterpolant(reshape(XN,[],1),reshape(YN,[],1),reshape(temp,[],1));
%         F = griddedInterpolant(XN,YN,temp,'linear');
        phiv = F(xv,yv); 
        
        phiv = naninterp(phiv,0);
    

end

end