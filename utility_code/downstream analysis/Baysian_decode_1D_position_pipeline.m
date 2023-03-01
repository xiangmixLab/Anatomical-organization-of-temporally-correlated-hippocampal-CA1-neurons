% Baysian decoding using codes from Guillaume et al. 2020, code aviliable
% at  https://github.com/etterguillaume/CaImDecoding.

function [decoded_position,ori_position_bin,decoded_accuracy,decoded_error]=Baysian_decode_1D_position_pipeline(neuron_train,behav_train,neuron_test,behav_test,binsize,cell_list)

%% pick cells
if ~isempty(cell_list)
    neuron_train=pickNeuronsByList(neuron_train,cell_list);
    neuron_test=pickNeuronsByList(neuron_test,cell_list);
end
%% Binarize calcium trace
sampling_frequency = 15; % This data set has been sampled at 15 images per second
% z_threshold = 2; % A 2 standard-deviation threshold is usually optimal to differentiate calcium ativity from background noise
binarized_trace_tr=[];
binarized_trace_ts=[];
for i=1:size(neuron_train.C,1)
    [binarized_trace_tr(:,i),~,ntrace] = extract_binary_nC(neuron_train.C(i,:), sampling_frequency, 0.1*max(neuron_train.C(i,:)./std(neuron_train.C(i,:))));
    [binarized_trace_ts(:,i)] = extract_binary_nC(neuron_test.C(i,:), sampling_frequency, 0.1*max(neuron_test.C(i,:)./std(neuron_test.C(i,:))));
end
%% Interpolate behavior
% In most cases, behavior data has to be interpolated to match neural temporal
% activity assuming it has a lower sampling frequency than neural activity
% interp_behav_vec_tr(:,1) = interpolate_behavior(behav_train.position(:,1), behav_train.time, neuron_train.time); % in the X dimension
% interp_behav_vec_tr(:,2) = interpolate_behavior(behav_train.position(:,2), behav_train.time, neuron_train.time); % in the Y dimension

[interp_behav_vec_tr,behav_tr_time]=align_neuron_behav(neuron_train,behav_train.time,behav_train.position);
[interp_behav_vec_ts,behav_ts_time]=align_neuron_behav(neuron_test,behav_test.time,behav_test.position);

%% behav to 1d
max_x=max(interp_behav_vec_tr(:,1))-min(interp_behav_vec_tr(:,1));
max_y=max(interp_behav_vec_tr(:,2))-min(interp_behav_vec_tr(:,2));
if max_x>max_y % hor track
    interp_behav_vec_tr(:,2)=ones(size(interp_behav_vec_tr,1),1);
    interp_behav_vec_ts(:,2)=ones(size(interp_behav_vec_ts,1),1);
else % vet track
    interp_behav_vec_tr(:,1)=ones(size(interp_behav_vec_tr,1),1);
    interp_behav_vec_ts(:,1)=ones(size(interp_behav_vec_ts,1),1);
end

%% Compute velocities
% In this example, we will ignore moments when the mouse is immobile
behav_struct_tr.position=interp_behav_vec_tr;
behav_struct_tr.time=behav_tr_time;
[ velocity_tr ] = behav_velo_cal(behav_struct_tr,[],'r');
behav_struct_ts.position=interp_behav_vec_ts;
behav_struct_ts.time=behav_ts_time;
[ velocity_ts ] = behav_velo_cal(behav_struct_ts,[],'r');

%% Find periods of immobility
% This will be usefull later to exclude immobility periods from analysis
min_speed_threshold = 5; % 5mm/sec
running_tr_ts = velocity_tr > min_speed_threshold;
running_ts_ts = velocity_ts > min_speed_threshold;

%% Compute occupancy and joint probabilities
% You can use 'min(behav_vec)' and 'max(behav_vec)'  to estimate the
% boundaries of the behavior vector (in our case, location in space)
bin_size = binsize;
X_bin_vector = min(interp_behav_vec_tr(:,1)):bin_size:max(interp_behav_vec_tr(:,1))+bin_size;
X_bin_centers_vector = X_bin_vector + bin_size/2;
X_bin_centers_vector(end) = [];

Y_bin_vector = min(interp_behav_vec_tr(:,2)):bin_size:max(interp_behav_vec_tr(:,2))+bin_size;
Y_bin_centers_vector = Y_bin_vector + bin_size/2;
Y_bin_centers_vector(end) = [];

%% next steps, decoding

decoded_position={};
decoded_accuracy=[];

prob_being_active=[];
tuning_map_data=[];
all_posterior={};
for cell_i = 1:size(binarized_trace_tr,2)
    [~,Posterior, occupancy_map, prob_being_active(cell_i), tuning_map_data(:,:,cell_i) ] = extract_2D_information(binarized_trace_tr(:,cell_i), interp_behav_vec_tr, X_bin_vector, Y_bin_vector, running_tr_ts);
    all_posterior{cell_i}=Posterior;
end

% Minimal a priori (use to remove experimental a priori)
[~,~, occupancy_map1,~,~] = extract_2D_information(binarized_trace_ts(:,cell_i), interp_behav_vec_ts, X_bin_vector, Y_bin_vector, running_ts_ts);
occupancy_vector = occupancy_map1./occupancy_map1*(1/length(occupancy_map1)); % actually, even using the occupancy_map in training is ok, no difference after minimize a priori

[decoded_probabilities] = bayesian_decode2D(binarized_trace_ts, occupancy_vector, prob_being_active, tuning_map_data, 1:size(binarized_trace_ts,2));
[filt_decoded_probabilities] = bayesian_temporal_filter2D_adapted(decoded_probabilities); % 5 sec smoothing

for p=1:size(filt_decoded_probabilities,3)
    dpt1=squeeze(filt_decoded_probabilities(:,:,p));
    dpt1(:,[1 size(dpt1,2)])=dpt1(:,[1 size(dpt1,2)])*0;
    dpt1([1 size(dpt1,1)],:)=dpt1([1 size(dpt1,1)],:)*0;
    filt_decoded_probabilities(:,:,p)=dpt1;
end

decoded_positiont=zeros(size(decoded_probabilities,3),2);
fpt_all={};
for p=1:size(decoded_probabilities,3)
    dpt1=squeeze(decoded_probabilities(:,:,p)-nanmean(decoded_probabilities,3));
%     fpt=filter2DMatrices_021521(dpt1,5,2);
%     fpt_all{p}=fpt{1}>0.7*max(fpt{1}(:));
%     fpt_all{p}=dpt1>0.8*max(dpt1(:));
%     stats=regionprops(fpt_all{p});
%     centroidd=[];
%     Area_all=[stats.Area];
%     major_area_idx=find(Area_all>0.5*max(Area_all));
%     
%     ctt=1;
%     centroidd=[];
%     for k=major_area_idx
%         centroidd(ctt,:)=stats(k).Centroid;
%         ctt=ctt+1;
%     end
%     
%     if ~isempty(centroidd)
%         decoded_positiont(p,:)=mean(centroidd,1);
%     else
%         decoded_positiont(p,:)=[nan nan];
%     end
    if ~isempty(dpt1)&&nansum(dpt1(:))~=0
        dpt_max=[];
        [dpt_max(:,2),dpt_max(:,1)]=find(dpt1==max(dpt1(:)));
        if p>1
            if ~isnan(sum(decoded_positiont(p-1,:)))
                distt=sum((dpt_max-decoded_positiont(p-1,:)).^2,2).^0.5;
                if size(distt,1)==1
                    [decoded_positiont(p,:)]=dpt_max(distt==min(distt),:);
                else
                    decoded_positiont(p,:)=nanmean(dpt_max,1);
                end
            else
                decoded_positiont(p,:)=nanmean(dpt_max,1);
            end
        else
            decoded_positiont(p,:)=nanmean(dpt_max,1);
        end
    else
        decoded_positiont(p,:)=[nan nan];
    end

end
decoded_positiont=fillmissing(decoded_positiont,'linear');
decoded_position=smoothdata(decoded_positiont,'SmoothingFactor',0.5);
ori_position_bin=smoothdata(round(interp_behav_vec_ts/binsize),'SmoothingFactor',0.05);


if max_x>max_y % hor track
    decoded_accuracy_t=corrcoef(ori_position_bin(:,1),decoded_position(:,1));
    decoded_accuracy=decoded_accuracy_t(2); % avg accuracy from x and y
    decoded_error=nanmean(abs(ori_position_bin(:,1)-decoded_position(:,1)));
else % vet track
    decoded_accuracy_t=corrcoef(ori_position_bin(:,2),decoded_position(:,2));
    decoded_accuracy=decoded_accuracy_t(2); % avg accuracy from x and y
    decoded_error=nanmean(abs(ori_position_bin(:,2)-decoded_position(:,2)));
end