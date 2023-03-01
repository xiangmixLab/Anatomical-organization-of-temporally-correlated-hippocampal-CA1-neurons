function quad_trace_layerAvg=quad_layer_mean(quad_trace,layerNum,nC_leng)

qvstack=[];
for i=1:layerNum
    q{i}=quad_trace(i:layerNum:end);
    if length(q{i})~=nC_leng
        if length(q{i})>nC_leng
            q{i}=q{i}(1:nC_leng);
        else
            q{i}=[q{i} zeros(nC_leng-length(q{i}))]; % should not be too much, maybe 1-2, except at the end of experiment, trim all quads after end of neuron
        end
    end
    qvstack=[qvstack;q{i}];    
end

quad_trace_layerAvg=nanmean(qvstack,1); % velo average across 3 layers
quad_trace_layerAvg=quad_trace_layerAvg>1;

