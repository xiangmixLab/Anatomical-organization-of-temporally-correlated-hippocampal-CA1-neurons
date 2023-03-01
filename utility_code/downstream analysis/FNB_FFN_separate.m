function [neuron_sep,behav_sep]=FNB_FFN_separate(neuron,behav,blocker_pos)

behavpos1=behav.position;

idx1=behavpos1(:,1)<=blocker_pos-5;
idx2=behavpos1(:,1)>=blocker_pos+5;

behav_sep{1}=behav;
behav_sep{2}=behav;
behav_sep{1}.position=behav_sep{1}.position(idx1,:);
behav_sep{2}.position=behav_sep{2}.position(idx2,:);
behav_sep{1}.time=behav_sep{1}.time(idx1,:);
behav_sep{2}.time=behav_sep{2}.time(idx2,:);

ntime = double(neuron.time);
ntime = ntime(1:2:end);
ntime = resample(ntime,size(neuron.C,2),length(ntime));
behavtime=behav.time;
idx_resample1=interp1(behavtime,double(idx1),ntime,'nearest');
idx_resample1_1=interp1(behavtime,double(idx1),neuron.time,'nearest');
idx_resample2=interp1(behavtime,double(idx2),ntime,'nearest');
idx_resample2_1=interp1(behavtime,double(idx2),neuron.time,'nearest');

idx_resample1(isnan(idx_resample1))=0;
idx_resample1_1(isnan(idx_resample1_1))=0;
idx_resample2(isnan(idx_resample2))=0;
idx_resample2_1(isnan(idx_resample2_1))=0;


idx_resample1=logical(idx_resample1);
idx_resample1_1=logical(idx_resample1_1);
idx_resample2=logical(idx_resample2);
idx_resample2_1=logical(idx_resample2_1);

neuron_sep{1}=neuron.copy;
neuron_sep{2}=neuron.copy;

if length(neuron.time)==size(neuron.C,2)
    neuron_sep{1}.time=neuron.time(idx_resample1_1);
else
    neuron_sep{1}.time=neuron.time(idx_resample1);
end
neuron_sep{1}.C=neuron.C(:,idx_resample1);
neuron_sep{1}.C_raw=neuron.C_raw(:,idx_resample1);
neuron_sep{1}.S=neuron.S(:,idx_resample1);
neuron_sep{1}.trace=[];
neuron_sep{1}.C_df=[];

if length(neuron.time)==size(neuron.C,2)
    neuron_sep{2}.time=neuron.time(idx_resample2_1);
else
    neuron_sep{2}.time=neuron.time(idx_resample2);
end
neuron_sep{2}.C=neuron.C(:,idx_resample2);
neuron_sep{2}.C_raw=neuron.C_raw(:,idx_resample2);
neuron_sep{2}.S=neuron.S(:,idx_resample2);
neuron_sep{2}.trace=[];
neuron_sep{2}.C_df=[];

