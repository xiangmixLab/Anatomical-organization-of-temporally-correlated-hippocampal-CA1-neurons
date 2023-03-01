function [group,CMf,Zf,cophList,cophList_adjusted]=cluster_determine_by_suoqin_NMF_firstPeakCoph(neuron0,N,maxKguess,optimalK_preDefine)
cophList=[];
cophList_adjusted=0;
tic;

neuron0.C(isnan(neuron0.C(:)))=0;
neuron0.S=C_to_peakS(neuron0.C);
        
if isempty(optimalK_preDefine)
    
    parfor k=2:maxKguess
        CM = consensusKmeans_adapted(neuron0,nanmean(2*nanstd(neuron0.S,[],2)),k,N,[],[]); % return the simimarity matrix between paired neurons
        [~,~,~,coph] = nmforderconsensus0(CM,k);
        cophList(k-1)=coph;
        disp(['got ',num2str(k)]);
    end
    toc;

    % cophList_adjusted= smoothdata(cophList,'SmoothingFactor',0.5);

    [~,ploc]=findpeaks(cophList);
    KList=[2:maxKguess];
    if ~isempty(ploc)
        optimalK=KList(min(ploc));
    else
        cophList_dif=diff(cophList);
        optimalK=KList(find(cophList_dif==min(cophList_dif)));
    end

else
    optimalK=optimalK_preDefine;
end
    
CMf = consensusKmeans_adapted(neuron0,nanmean(2*nanstd(neuron0.S,[],2)),optimalK,N,[],[]);
Zf = linkage(CMf,'complete');
group=cluster(Zf,'maxclust',optimalK);

% remove neurons with overall very small correlation with other neurons
corrAll=squareform(1-pdist(neuron0.C,'correlation'));
maxCorrAll=max(corrAll);
idx_smallCorr=maxCorrAll<0.3