function [phiv] = interpNCOM_grid_parallel(XN,YN,phiN,xv,yv)
% this function takes NCOM xN,yN,zN,phi data
% interpolates in xy using parfor to obtain phi (T,S,U,V)
% this function assumes variables are on the same z grid

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
        indx = ~isnan(temp);
%         temp = reshape(temp,[],1);
        if sum(sum(indx))>2
                F = scatteredInterpolant(reshape(XN(indx),[],1),reshape(YN(indx),[],1),...
                    reshape(temp(indx),[],1),'linear','nearest');
%                 F = griddedInterpolant(XN,YN,squeeze(phiN(:,:,i)),'linear');
                phiv(:,i) = F(xv,yv); 
        else
            phiv(:,i) = xv+nan;

        end
    end
    
    % now interp in z, extrapolate unnknown values
%     phiv = interp1(zN,phiv_z',zv,'linear','extrap')';
    
    phiv = naninterp(phiv,0);
    
elseif  length(size(phiN))==2 % 2d variable
       temp = squeeze(phiN);
%         temp = inpaintn(squeeze(phiN));
        
        F = griddedInterpolant(XN,YN,temp,'linear');
        phiv = F(xv,yv); 
        
        phiv = naninterp(phiv,0);
    

end

end