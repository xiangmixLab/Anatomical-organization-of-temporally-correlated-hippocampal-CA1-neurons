function [MI, posterior, occupancy_map, prob_being_active] = extract_2D_information_adapted(binarized_trace, count, countTime_pt)
%MSPLACE_CELL_OF_POSITIONS Analyse and plot spatial information
%   This function plots in space the position related to calcium

% Copyright (C) 2017-2019 by Guillaume Etter
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or any
% later version.
% Contact: etterguillaume@gmail.com

% binary_trace: logical vector representing neurons active/inactive periods
%
% interp_behav_vec: m x 2  matrix representing behavioral state in two dimensions (eg position in a
% maze). Has to be the same size as binary_trace
%
% ca_time: time vector for both binary_trace and interp_behav_trace. If time correspondance
% is different between these two variables, you can use the fonction
% interp_behav to interpolate the behavioral trace
%
% OPTIONAL vectors:
% inclusion_vec: logical vector including (1) or excluding (0) corresponding timestamps. Has to be the same size as binary_trace.

%% Ignore excluded periods

%% Create bin vectors

%% Compute joint probabilities (of cell being active while being in a state bin)
MI = 0;

prob_being_active=sum(binarized_trace)/length(binarized_trace);
occupancy_map=countTime_pt./length(binarized_trace);

activity_in_bin=count;
inactivity_in_bin=countTime_pt-count;

likelihood=activity_in_bin./countTime_pt;

MI_list=[];
for y = 1:size(countTime_pt,1)
    for x = 1:size(countTime_pt,2)     
        
        if countTime_pt(y,x)~=0
            joint_prob_active = activity_in_bin(y,x)/length(binarized_trace);% count/total time
            joint_prob_inactive = inactivity_in_bin(y,x)/length(binarized_trace);
            prob_in_bin = occupancy_map(y,x);

            if joint_prob_active ~= 0
                MI = MI + joint_prob_active*log2(joint_prob_active./(prob_in_bin*prob_being_active));
                MI_list=[MI_list,joint_prob_active*log2(joint_prob_active./(prob_in_bin*prob_being_active))];
            end
            if joint_prob_inactive ~= 0
                MI = MI + joint_prob_inactive*log2(joint_prob_inactive./(prob_in_bin*(1-prob_being_active)));
                MI_list=[MI_list,joint_prob_inactive*log2(joint_prob_inactive./(prob_in_bin*(1-prob_being_active)))];
            end
        end
    end
end

posterior = likelihood.*occupancy_map/prob_being_active;

end



