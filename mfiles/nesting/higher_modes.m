function [w,c] = higher_modes(z,N,omega,nmode,topbc)

     g = 9.81;
     if(diff(z)<0)
       error('z must be monotonically increasing');
     end

     Nk = length(N);
     dz = zeros(Nk,1);
     dziph = zeros(Nk+1,1);

     dziph(2:end-1) = z(2:end)-z(1:end-1);
     dziph(1) = 2*dziph(2)-dziph(3);
     dziph(end) = 2*dziph(end-1)-dziph(end-2);

     % a(k)*w(k-1)+b(k)*w(k)+c(k)*w(k+1) + f/c^2*w(k) = 0     
     a = zeros(Nk,1);
     b = zeros(Nk,1);
     c = zeros(Nk,1);
    
     dz = 0.5*(dziph(1:end-1)+dziph(2:end));
     a = 1./(dziph(1:end-1).*dz);
     b = -(1./dziph(1:end-1)+1./dziph(2:end))./dz;
     c = 1./(dziph(2:end).*dz);

     f = (N.^2-omega^2);
     a = -a./f;     
     b = -b./f;
     c = -c./f;

     % Bottom no-flux w=0.
     b(1) = b(1)-a(1);

     if(strcmp(topbc,'rigid lid'))
       % Top rigid lid: w=0 
       b(Nk) = b(Nk)-c(Nk);
     elseif(strcmp(topbc,'free surface'))
       % Free-surface: dw/dz + g/c^2 w = 0
       a(Nk) = -1/(g*dziph(end-1));
       b(Nk) = 1/(g*dziph(end-1));
     elseif(strcmp(topbc,'free surface kundu'))
       % Free-surface: dw/dz + N^2/g w = 0
       alpha =  (2 - dziph(Nk+1)*(N(Nk)^2-omega^2)/g)/...
                (2 + dziph(Nk+1)*(N(Nk)^2-omega^2)/g);
       b(Nk) = b(Nk)+alpha*c(Nk);
     else
       error('topbc input must be either ''rigid lid'' or ''free surface''');
     end

     A = diag(a(2:end),-1)+diag(b,0)+diag(c(1:end-1),1);
     [X,L] = eig(A,'nobalance');

     [lams,is] = sort(diag(L));
     X = X(:,is);

     w = X(:,[1:nmode]);
     c = sqrt(1./lams([1:nmode]));



