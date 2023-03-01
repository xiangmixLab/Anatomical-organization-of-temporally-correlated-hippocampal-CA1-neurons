function [nnew,behavpos_new,behavtime_new,idx,idx_resample,idx_resample1,velocity_new,velocity_resample1]=linearTrack_remove_water(neuron,behavpos,behavtime,watersize,smoothingFactor)

%% get the maximal position of x and y axis, to determine how the track is oriented
max_behavpos_x=max(behavpos(:,1));
max_behavpos_y=max(behavpos(:,2));

try
    nnew=neuron.copy;
catch
    nnew=Sources2D;
    nnew.A=neuron.A;
    nnew.C=neuron.C;
    nnew.C_raw=neuron.C_raw;
    nnew.S=C_to_peakS(nnew.C);
    nnew.time=neuron.time; % simple version
    nnew.centroid=neuron_centroid_calculation(neuron,[240,376]);
    
    try
        nnew.Cn=neuron.Cn;
    catch
    end
    try
        nnew.imageSize=neuron.imageSize;
    catch
        nnew.imageSize=[240,376];
    end
end

% horizontal track
if max_behavpos_x>max_behavpos_y
    % define watersize% of the track length as water size
    water_size=(max(behavpos(:,1))-min(behavpos(:,1)))*watersize;
    % water area is defined as end+/-water size
    min_pos=min(behavpos(:,1))+water_size;
    max_pos=max(behavpos(:,1))-water_size;
    
    % get the position index that between the two water region
    behavpos1=behavpos;
    diff_behavpos1=max(behavpos1,[],1)-min(behavpos,[],1);
    behavpos1(:,diff_behavpos1==min(diff_behavpos1))=1;
    velocity=extract_velocity_adapted(behavpos1, behavtime/1000,smoothingFactor);
    idx=logical((behavpos(:,1)>=min_pos).*(behavpos(:,1)<=max_pos));
    behavpos_new=behavpos(idx,:);
    behavtime_new=behavtime(idx);
    velocity_new=velocity(idx);
    
    ntime = double(nnew.time);
    ntime = ntime(1:2:end);
    ntime = resample(ntime,size(nnew.C,2),length(ntime));
    
    idx_resample=interp1(behavtime,double(idx),ntime,'nearest');
    idx_resample1=interp1(behavtime,double(idx),nnew.time,'nearest');
    
    idx_resample(isnan(idx_resample))=0;
    idx_resample1(isnan(idx_resample1))=0;
    
    idx_resample=logical(idx_resample);
    idx_resample1=logical(idx_resample1);
    
    nnew.C=nnew.C(:,idx_resample);
    nnew.C_raw=nnew.C_raw(:,idx_resample);
    nnew.S=nnew.S(:,idx_resample);
    nnew.time=nnew.time(idx_resample1);
    
    velocity_resample=interp1(behavtime,velocity,ntime,'nearest');
    velocity_resample1=velocity_resample(idx_resample);
else
    %vertical track
    water_size=(max(behavpos(:,2))-min(behavpos(:,2)))*0.1;
    min_pos=min(behavpos(:,2))+water_size;
    max_pos=max(behavpos(:,2))-water_size;
    
    behavpos1=behavpos;
    diff_behavpos1=max(behavpos1,[],1)-min(behavpos,[],1);
    behavpos1(:,diff_behavpos1==min(diff_behavpos1))=1;

    velocity=extract_velocity_adapted(behavpos1, behavtime/1000,smoothingFactor);
    idx=logical((behavpos(:,2)>=min_pos).*(behavpos(:,2)<=max_pos));
    behavpos_new=behavpos(idx,:);
    behavtime_new=behavtime(idx);
    velocity_new=velocity(idx);
    
    ntime = double(nnew.time);
    ntime = ntime(1:2:end);
    ntime = resample(ntime,size(nnew.C,2),length(ntime));
    
    idx_resample=interp1(behavtime,double(idx),ntime,'nearest');
    idx_resample1=interp1(behavtime,double(idx),nnew.time,'nearest');
    
    idx_resample(isnan(idx_resample))=0;
    idx_resample1(isnan(idx_resample1))=0;
    
    idx_resample=logical(idx_resample);
    idx_resample1=logical(idx_resample1);  
    
    nnew.C=nnew.C(:,idx_resample);
    nnew.C_raw=nnew.C_raw(:,idx_resample);
    nnew.S=nnew.S(:,idx_resample);
    nnew.time=nnew.time(idx_resample1);
    
    velocity_resample=interp1(behavtime,velocity,ntime,'nearest');
    velocity_resample1=velocity_resample(idx_resample);
end