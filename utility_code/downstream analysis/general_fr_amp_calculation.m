function [fr_all,amp_all,fr_t,amp_t]=general_fr_amp_calculation(neuronIndividuals_new,del_idx,low_fr_thresh)
% the conditions of neuronIndividuals_new has to be done prior
% neuron group has to be done prior
%% merge from fig3: fr/amp per cell
fr_k={};
amp_k={};

for i=1:length(neuronIndividuals_new)
    
    nC=neuronIndividuals_new{i}.C;     
    nC(del_idx,:)=[];
    

    % to dF/F
    for j=1:size(nC,1)
        t=nC(j,:);
        nC(j,:)=(t-nanmean(t(t>0)))./nanmean(t(t>0)); % long slience period may exist
    end
    thresh=0.1*nanmax(nC,[],2);
    
    if max(nC(:))>0
        %get peak amplitude from dF/F nC
        nS=C_to_peakS(nC);
        
        % threshold peak
        for k=1:size(nS,1)
            St=nS(k,:);
            St(St<thresh(k,1))=0;
            nS(k,:)=St;
        end

        fr=nansum(nS>0,2)./(size(nS,2)/15);
        amp=nansum(nS,2)./nansum(nS>0,2);
        fr(isnan(fr))=0;
        amp(isnan(amp))=0;

        % thres fr
        if ~isempty(low_fr_thresh)
            low_fr=fr<low_fr_thresh;
            fr(low_fr)=0;
            amp(low_fr)=0;
        end

        fr_k{i}=fr;
        amp_k{i}=amp;
    else
        fr_k{i}=zeros(size(nC,1),1);
        amp_k{i}=zeros(size(nC,1),1);
    end
end

fr_t=cell2mat(fr_k);
amp_t=cell2mat(amp_k);

larger_than_0_fr_cir=sum(fr_t>0,2);
larger_than_0_amp_cir=sum(amp_t>0,2);    

fr_all=sum(fr_t,2)./(larger_than_0_fr_cir);
amp_all=sum(amp_t,2)./(larger_than_0_amp_cir);