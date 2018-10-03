function [PROF] = get_profile(dirname)
%
% Plot profiles of u(z,t) at the different stations specified in
% ../rundata/dataxy.dat
%
% dirname = '../data';
EMPTY = 999999;

fname = [dirname,'/profdata.dat'];


fid = fopen(fname,'rb');
if fid==-1 % there is no file
    PROF=[];
    return
else
numTotalDataPoints = fread(fid,1,'int32');
numInterpPoints = fread(fid,1,'int32');
Nkmax = fread(fid,1,'int32');
nsteps = fread(fid,1,'int32');
ntoutProfs = fread(fid,1,'int32');  
dt = fread(fid,1,'float64');
dz = fread(fid,Nkmax,'float64');
dataIndices = fread(fid,numTotalDataPoints,'int32');
dataXY = fread(fid,2*numTotalDataPoints,'float64');
xv = reshape(fread(fid,numInterpPoints*numTotalDataPoints,'float64'),numInterpPoints,numTotalDataPoints);
yv = reshape(fread(fid,numInterpPoints*numTotalDataPoints,'float64'),numInterpPoints,numTotalDataPoints);
fclose(fid);

starttime = getvalue([dirname,'/suntans.dat'],'starttime');

% load profile names
dataNames = dir([dirname, '/dataxy_names.dat']);
if ~isempty(dataNames) % file exists   
dataNames = importdata([dirname, '/dataxy_names.dat']);
end


%%
dataX = dataXY(1:2:end);
dataY = dataXY(2:2:end);

z = getz(dz);

fname = [dirname,'/u.data.prof'];
fid = fopen(fname,'rb');
if fid==-1 % there is no file
    PROF=[];
    return
else

    data = fread(fid,'float64');
    data(find(data==EMPTY))=nan;
    fclose(fid);
    
if isempty(data)
    PROF=[];
    return
end


    nout = length(data)/(3*Nkmax*numInterpPoints*numTotalDataPoints);
    udata = reshape(data,Nkmax,numInterpPoints,numTotalDataPoints,3,nout);
    udata = squeeze(udata);

    Time = ones(Nkmax,1)*[1:nout]*dt*ntoutProfs; % s since start
    mtime = Time/24/3600+datenum(sprintf('%.6f',starttime),'yyyymmdd.HHMMSS');

    Z = z*ones(1,nout);
    Tday = 86400;

    u = squeeze(udata(:,:,1,:));
    v = squeeze(udata(:,:,2,:));
    w = squeeze(udata(:,:,3,:));

    for loc=1:size(u,2)
      uplot{loc} = squeeze(u(:,loc,:));
      vplot{loc} = squeeze(v(:,loc,:));
      wplot{loc} = squeeze(w(:,loc,:));
      kmax = max(find(~isnan(uplot{loc}(:,1))));
      d0(loc) = sum(dz(1:kmax));
      u0{loc} = sum((dz(1:kmax)*ones(1,nout)).*uplot{loc}(1:kmax,:))/d0(loc);
      v0{loc} = sum((dz(1:kmax)*ones(1,nout)).*vplot{loc}(1:kmax,:))/d0(loc);
      w0{loc} = sum((dz(1:kmax)*ones(1,nout)).*wplot{loc}(1:kmax,:))/d0(loc);
    end

    PROF.uplot = uplot;
    PROF.vplot = vplot;
    PROF.wplot = wplot;
    PROF.u0 = u0;
    PROF.v0 = v0;
    PROF.w0 = w0;
    PROF.d0 = d0;
    PROF.dz = dz;
    PROF.Time = Time;
    PROF.mtime = mtime;
    PROF.Tday = Tday;
    PROF.z = z;
    PROF.Z = Z;
    PROF.dataX = dataX;
    PROF.dataY = dataY;
    PROF.dataNames=dataNames;

    %% other variables

    % temperature
    fname = [dirname,'/T.data.prof'];
    fid = fopen(fname,'rb');
    if fid==-1 % there is no file
        PROF.Tplot=[];   
    else
        data = fread(fid,'float64');
        data(find(data==EMPTY))=nan;
        fclose(fid);
        Tdata = reshape(data,Nkmax,numInterpPoints,numTotalDataPoints,nout);
        Tdata = squeeze(Tdata);
        for loc=1:size(Tdata,2)
          PROF.Tplot{loc} = squeeze(Tdata(:,loc,:));
        end
    end

    % salinity
    fname = [dirname,'/s.data.prof'];
    fid = fopen(fname,'rb');
    if fid==-1 % there is no file
        PROF.Splot=[];   
    else
        data = fread(fid,'float64');
        data(find(data==EMPTY))=nan;
        fclose(fid);
        Sdata = reshape(data,Nkmax,numInterpPoints,numTotalDataPoints,nout);
        Sdata = squeeze(Sdata);
        for loc=1:size(Sdata,2)
          PROF.Splot{loc} = squeeze(Sdata(:,loc,:));
        end
    end

    % nh pressure
    fname = [dirname,'/q.data.prof'];
    fid = fopen(fname,'rb');
    if fid==-1 % there is no file
        PROF.qplot=[];   
    else
        data = fread(fid,'float64');
        data(find(data==EMPTY))=nan;
        fclose(fid);
        qdata = reshape(data,Nkmax,numInterpPoints,numTotalDataPoints,nout);
        qdata = squeeze(qdata);
        for loc=1:size(qdata,2)
          PROF.qplot{loc} = squeeze(qdata(:,loc,:));
        end
    end

    % free surface
    fname = [dirname,'/fs.data.prof'];
    fid = fopen(fname,'rb');
    if fid==-1 % there is no file
        PROF.etaplot=[];   
    else
        data = fread(fid,'float64');
        data(find(data==EMPTY))=nan;
        fclose(fid);
        etadata = reshape(data,numInterpPoints,numTotalDataPoints,nout);
        etadata = squeeze(etadata);
        for loc=1:size(etadata,1)
          PROF.etaplot{loc} = squeeze(etadata(loc,:));
        end
    end

    end
    %%
    % figure(1);
    % for loc=1:3
    %   
    %   uplot = squeeze(u(:,loc,:));
    %   kmax = max(find(~isnan(uplot(:,1))));
    %   d0 = sum(dz(1:kmax));
    %   u0 = sum((dz(1:kmax)*ones(1,nout)).*uplot(1:kmax,:))/d0;
    %   
    %   subplot(3,1,loc)
    %   pcolor(Time/Tday,Z,uplot-ones(Nkmax,1)*u0);
    %   title(sprintf('Baroclinic u(z,t) at Location %d (x = %.0f km)',loc,dataX(loc)/1000));
    %   shading flat;
    %   axis([0 max(Time(1,:)/Tday) -d0 0]);
    %   colorbar;
    %   
    %   if(loc==3) 
    %     xlabel('Time (days)'); 
    %   else
    %     set(gca,'xticklabel','');
    %   end
    %   ylabel('Depth (m)');
    % end
end