function [gp2,gp2_highcorr]=cluster_idx_matching_ensembleRatemap(group1,group2,n1,n2,b1,b2)

[ratemap1,~,~,countTime1] = calculatingCellSpatialForSingleData_Suoqin(n1,b1.position,b1.time,b1.ROI,10,1:size(n1.S,1),3*std(n1.S,[],2),'S',[],[],[0.1 1000000],10);
[ratemap2,~,~,countTime2] = calculatingCellSpatialForSingleData_Suoqin(n2,b2.position,b2.time,b2.ROI,10,1:size(n2.S,1),3*std(n2.S,[],2),'S',[],[],[0.1 1000000],10);

[group12,group21]=determineReorganizedCluster_072922(group1,group2);
[esm_ratemap_all]=determineReorganizedEnsembleRatemap_072922(group1,group2,group12,group21,ratemap1,ratemap2);

ensembleMap1=esm_ratemap_all(1,:);
ensembleMap2=esm_ratemap_all(2,:);

corrMat=[];
for i=1:length(ensembleMap1)
    for j=1:length(ensembleMap2)
        map1=ensembleMap1{i};
        map2=ensembleMap2{j};
        
        corr=rateMap_correlation({map1},{map2},countTime1,countTime2,1,0);
        
        corrMat(i,j)=corr;
    end
end

% assign idx
correspondence=[];

% 1. assign correspondence if c1 from cluster 1, c2 from cluster 2, c1 have
% max corr with c2 and c2 has max corr with c1 compared to other clusters
for i=1:size(corrMat,1)
    if max(corrMat(i,:))>0 % handle non equal number of clusters
        max_currRow=find(corrMat(i,:)==max(corrMat(i,:)));
        corresponding_column=corrMat(:,max_currRow);
        if find(corresponding_column==max(corresponding_column))==i && corrMat(i,max_currRow)>=0.6
            correspondence=[correspondence;i,max_currRow];
        end
    end
end

% 3. assign highcorr correspondence
gp2_highcorr=group2*0;
for i=1:size(correspondence,1)
    gp2_highcorr(group2==correspondence(i,2))=correspondence(i,1);
end
gp2_highcorr(gp2_highcorr==0)=-1;

% 2. for remaining ones
ugp=unique(group1);
while size(correspondence,1)<length(ugp)
    for i=1:size(corrMat,1)
        if ~ismember(i,correspondence(:,1)) % if already assigned, jump; if not , do
            [~,idx]=sort(corrMat(i,:),'ascend');
            for j=1:length(idx)
                if ismember(idx(j),correspondence(:,2))
                    continue;
                end
                % reach current unused maximum corr
                corresponding_column=corrMat(:,idx(j));
                corresponding_column(correspondence(:,1))=-1;
                containment=find(corresponding_column==max(corresponding_column));
                if ismember(i,containment)
                    correspondence=[correspondence;i,idx(j)];
                end
            end
        end
    end
end

% 3. assign correspondence
gp2=group2*0;
for i=1:size(correspondence,1)
    gp2(group2==correspondence(i,2))=correspondence(i,1);
end
gp2(gp2==0)=-1;
% test
% illustrating_ensemble_cluster_ratemap_suppl25({{{n1,n2}}},{{{b1,b2}}},[1,1],[1,1],[1,2],{{group1,gp2}})
