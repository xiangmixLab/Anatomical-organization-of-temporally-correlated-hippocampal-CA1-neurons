function [binarized_trace, filtered_trace,norm_trace,d1_trace] = extract_binary_nC(calcium_trace,sampling_frequency, z_threshold)
%extract_binary Converts raw calcium traces into binary traces
%   This function converts raw calcium traces inputs into a binarized output vector

% Copyright (C) 2017-2019 by Guillaume Etter
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or any
% later version.
% Contact: etterguillaume@gmail.com

% calcium_trace: double/float vector representing activity of a single neuron 
% 
% sampling frequency: frequency at which calcium imaging has beer fps)
%
% z_threshold: standard deviation threshold above which calcium activity is
% considered higher than noise

%% Parameters
% [bFilt,aFilt] = butter(2,  2/(sampling_frequency/2), 'low');

filtered_trace = calcium_trace;
norm_trace = filtered_trace./std(filtered_trace);
% norm_trace=norm_trace/max(norm_trace(:)); % comment out 061521 to use
% 2std threshold
d1_trace = diff(filtered_trace);
d1_trace(end+1) = 0;

binarized_trace = calcium_trace*0;
% binarized_trace(norm_trace>z_threshold & d1_trace>0) = 1;
binarized_trace(norm_trace>z_threshold) = 1; % both side
end

