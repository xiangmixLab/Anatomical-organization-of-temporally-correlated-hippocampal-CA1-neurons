function [neuron_train,neuron_test,behav_train,behav_test]=decoding_split_trial_cell(neuron,behav,train_idx)

neuron_train.C=[];
neuron_train.C_raw=[];
neuron_train.S=[];
neuron_train.time=[0];

behav_train.position=[];
behav_train.time=[0];
behav_train.ROI=[];

neuron_test.C=[];
neuron_test.C_raw=[];
neuron_test.S=[];
neuron_test.time=[0];

behav_test.position=[];
behav_test.time=[0];
behav_test.ROI=[];
    
if ~isempty(neuron)
    neuron_train_cell=neuron(train_idx);
    neuron_test_cell=neuron(setdiff(1:length(neuron),train_idx));
    
    for i=1:length(neuron_train_cell)
        neuron_train.C=[neuron_train.C,neuron_train_cell{i}.C];
        neuron_train.C_raw=[neuron_train.C_raw,neuron_train_cell{i}.C_raw];
        neuron_train.time=[neuron_train.time;neuron_train_cell{i}.time+neuron_train.time(end)+nanmean(diff(neuron_train_cell{i}.time))];
    end
    neuron_train.time=neuron_train.time(2:end);
    neuron_train.S=C_to_peakS(neuron_train.C);
    
    
    for i=1:length(neuron_test_cell)
        neuron_test.C=[neuron_test.C,neuron_test_cell{i}.C];
        neuron_test.C_raw=[neuron_test.C_raw,neuron_test_cell{i}.C_raw];
        neuron_test.time=[neuron_test.time;neuron_test_cell{i}.time+neuron_test.time(end)+nanmean(diff(neuron_test_cell{i}.time))];
    end
    neuron_test.time=neuron_test.time(2:end);
    neuron_test.S=C_to_peakS(neuron_test.C);
end

if ~isempty(behav)
    behav_train_cell=behav(train_idx);
    behav_test_cell=behav(setdiff(1:length(behav),train_idx));
    
    behavROItr_stack=[];
    for i=1:length(neuron_train_cell)
        behav_train.position=[behav_train.position;behav_train_cell{i}.position];
        behav_train.time=[behav_train.time;behav_train_cell{i}.time+behav_train.time(end)+nanmean(diff(behav_train_cell{i}.time))];
        behavROItr_stack=[behavROItr_stack;behav_train_cell{i}.ROI];
    end
    behav_train.time=behav_train.time(2:end);
    behav_train.ROI=nanmean(behavROItr_stack,1);
    
    
    behavROIts_stack=[];
    for i=1:length(neuron_train_cell)
        behav_test.position=[behav_test.position;behav_test_cell{i}.position];
        behav_test.time=[behav_test.time;behav_test_cell{i}.time+behav_test.time(end)+nanmean(diff(behav_test_cell{i}.time))];
        behavROIts_stack=[behavROIts_stack;behav_test_cell{i}.ROI];
    end
    behav_test.time=behav_test.time(2:end);
    behav_test.ROI=nanmean(behavROIts_stack,1);
end
    
