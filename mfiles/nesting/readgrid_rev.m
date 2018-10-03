function [i,repeats,xall,yall,depthall,z,triall,Nc,Nkall] = readgrid_rev(dirname_g,numprocs)
% dirname_g='../data';
% numprocs=16;

  z = getz(load([dirname_g,'/vertspace.dat']));
   
  cp = load([dirname_g,'/cells.dat']);
  triall = cp(:,3:5);
  
  nn=1;
  for n=0:numprocs-1
    
    cdp = load([dirname_g,'/celldata.dat.',num2str(n)]);
   
    d{nn} = cdp(:,5);
    xv{nn} = cdp(:,2);
    yv{nn} = cdp(:,3);
    Nc(nn) = length(d{nn});
    Nk{nn} = cdp(:,6);
    
    nn=nn+1;
  end

  for n=0:numprocs-1
    if(n==0)
      xall = xv{1};
      yall = yv{1};
      depthall = d{1};
%       triall = tri{1};
      Nkall = Nk{1};
    else
      xall(end+1:end+Nc(n+1)) = xv{n+1};
      yall(end+1:end+Nc(n+1)) = yv{n+1};
      depthall(end+1:end+Nc(n+1)) = d{n+1};
%       triall(end+1:end+Nc(n+1),:) = tri{n+1};
      Nkall(end+1:end+Nc(n+1)) = Nk{n+1};
    end
  end

  R = sqrt(xall.^2+yall.^2);
  [dR,i] = sort(R);
  dR = dR(2:end)-dR(1:end-1);
  repeats = find(dR==0);
  xall(i(repeats)) = [];
  yall(i(repeats)) = [];
  depthall(i(repeats)) = [];
  triall(i(repeats),:) = [];
  Nkall(i(repeats)) = [];
