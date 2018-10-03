function [phiv] = interpNCOM(XN,YN,zN,phiN,xv,yv,zv)
% this function takes NCOM xN,yN,zN,phi data
% interpolates in xy to obtain phi (T,S,U,V)

% XN = NC.X;
% YN = NC.Y;
% zN = NC.z;
% phiN = NC.temp;
% zv=z_r;

XN = reshape(XN,[],1);
YN = reshape(YN,[],1);



if length(size(phiN))==3 %3d variable
    
    % interp in xy
    for i=1:size(phiN,3)
        temp = squeeze(phiN(:,:,i));
        temp = reshape(temp,[],1);
        
        F = scatteredInterpolant(XN,YN,temp);
        phiv_z(:,i) = F(xv,yv);   
    end
    
    % now interp in z
    phiv = interp1(zN,phiv_z',zv)';
    
    
elseif  length(size(phiN))==2 % 2d variable
    
        temp = reshape(phiN,[],1);
        
        F = scatteredInterpolant(XN,YN,temp);
        phiv = F(xv,yv); 
    

end

end