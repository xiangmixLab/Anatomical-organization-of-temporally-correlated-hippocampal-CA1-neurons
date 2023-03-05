% k-mean based consensus clustering

% INPUT: 
% neuron0: neuron data contains at least C (denoised calcium), and S (deconvoluted spikes)
% N: times for k-mean bootstraping
% maxKguess: upperbound for cluster number K if use auto cluster number
% search
% optimalK_preDefine: manually defined cluster number for the cluster
% calculation. If defined, auto cluster number search will be disabled

% OUTPUT:
% group: cluster result with pairs have correlation <= 0.3 excluded
% gp_ori: cluster result with all neuron included (used in paper)

% 2019-2023, Lujia Chen and Suoqin Jin, UC Irvine

function [group,gp_ori]=cluster_determine_by_suoqin_NMF_firstPeakCoph_022422(neuron0,N,maxKguess,optimalK_preDefine)
cophList=[];
cophList_adjusted=0;


neuron0.C(isnan(neuron0.C(:)))=0;
amp=C_to_peakS(neuron0.C);
neuron0.C(neuron0.C<3*std(amp,[],2))=0; % thresholding

if isempty(optimalK_preDefine)
    
    for k=2:maxKguess
        CM = consensusKmeans_adapted_latest(neuron0,k,N); % return the simimarity matrix between paired neurons
        [~,~,~,coph] = nmforderconsensus0(CM,k);
        cophList(k-1)=coph;
        disp(['got ',num2str(k)]);
    end

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
    
CMf = consensusKmeans_adapted_latest(neuron0,optimalK,N);
Zf = linkage(CMf,'complete');
group=cluster(Zf,'maxclust',optimalK);
gp_ori=group;
% remove neurons with overall very small correlation with other neurons
corrAll=squareform(1-pdist(neuron0.C,'correlation'));
maxCorrAll=max(corrAll);
idx_smallCorr=maxCorrAll<=0.3;
group(idx_smallCorr)=-1;
