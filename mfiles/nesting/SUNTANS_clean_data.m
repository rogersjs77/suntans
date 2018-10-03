% function finds ./data directories and deletes file

clear

% look in existing folder
d2 = dir('./data*');

if ~isempty(d2)
    
    rmdir('./data','s')
    disp(['removed ' d2.name '/data'])
end

d1 = dir('./');

for i=3:length(d1)
   
    % cycle through folders
    if d1(i).isdir
        cd([d1(i).folder '/' d1(i).name])
        delete('./mfiles/SUNTANS_results.mat');
%         delete('./mfiles/SUNTANS_quad_results.mat');
%         
        d2 = dir('./data*');
        if ~isempty(d2)
            rmdir('./data','s')
            disp(['removed ' d1(i).name '/data'])
        end
        cd '../'   
   
    end    
    
end
clear d1 d2