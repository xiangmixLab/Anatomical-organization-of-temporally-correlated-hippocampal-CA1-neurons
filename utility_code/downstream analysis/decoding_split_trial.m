function [neuron_train,neuron_test,behav_train,behav_test]=decoding_split_trial(neuron,behav,train_fraction)

if isempty(train_fraction)
    train_fraction=0.5;
end

neuron_train.C=[];
neuron_train.C_raw=[];
neuron_train.S=[];
neuron_train.time=[];

behav_train.position=[];
behav_train.time=[];
behav_train.ROI=[];

neuron_test.C=[];
neuron_test.C_raw=[];
neuron_test.S=[];
neuron_test.time=[];

behav_test.position=[];
behav_test.time=[];
behav_test.ROI=[];
    
if ~isempty(neuron)
    neuron_train.A=neuron.A;
    neuron_train.C=neuron.C(:,1:round(size(neuron.C,2)*train_fraction));
    neuron_train.C_raw=neuron.C_raw(:,1:round(size(neuron.C_raw,2)*train_fraction));
    neuron_train.S=C_to_peakS(neuron_train.C);
    neuron_train.time=neuron.time(1:round(length(neuron.time)*train_fraction));
    
    neuron_test.A=neuron.A;
    neuron_test.C=neuron.C(:,round(size(neuron.C,2)*train_fraction)+1:end);
    neuron_test.C_raw=neuron.C_raw(:,round(size(neuron.C,2)*train_fraction)+1:end);
    neuron_test.S=C_to_peakS(neuron_test.C);
    neuron_test.time=neuron.time(round(length(neuron.time)*train_fraction)+1:end)-neuron.time(round(length(neuron.time)*train_fraction)+1);
end

if ~isempty(behav)
    behav_train.position=behav.position(1:round(size(behav.position,1)*train_fraction),:);
    behav_train.time=behav.time(1:round(size(behav.time,1)*train_fraction));
    behav_train.ROI=behav.ROI;
    
    behav_test.position=behav.position(round(size(behav.position,1)*train_fraction)+1:end,:);
    behav_test.time=behav.time(round(size(behav.time,1)*train_fraction)+1:end)-behav.time(round(size(behav.time,1)*train_fraction)+1);
    behav_test.ROI=behav.ROI;
end
    
