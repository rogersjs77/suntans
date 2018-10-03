% a_SUNTANS_process
% Justin Rogers
% Stanford
clear, close all 

iprocess=1;
iplot = 1;
ireduce = 1;
quad=0; % process quad file
istep=1; % data points to skip
Tavg = 4*12.42*3600; % time from end of file to avg
datadir='../data'; % data directory

%%
save temp
if iprocess
    SUNTANS_process(datadir,istep,Tavg,quad)   
end

load temp
if ireduce
   a_SUNTANS_reduce 
   a_SUNTANS_reduce_3Dmovie
end

load temp
if iplot
    if quad
       SUNTANS_plot_quad 
    else
       SUNTANS_plot
    end
end



delete temp.mat
    