function [filt_decoded_probabilities] = bayesian_temporal_filter2D_adapted(decoded_probabilities)
%BAYESIAN_TEMPORAL_FILTER Summary of this function goes here
%   Detailed explanation goes here
Fs=15;
filt_decoded_probabilities=decoded_probabilities;
for i=1:size(decoded_probabilities,1)
    for j=1:size(decoded_probabilities,2)
        filt_decoded_probabilities(i,j,:)=smoothdata(squeeze(decoded_probabilities(i,j,:)),'movmean',Fs*5);
    end
end

