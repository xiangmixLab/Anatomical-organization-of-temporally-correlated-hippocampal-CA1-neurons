function [all_crossCorr_intra,all_crossCorr_intra_border,all_crossCorr_inter_border]=cluster_crossCorr_border(group,foldername,maxLagg)

% 1. calculating boundary idx
for i=1:size(group,1)
    load([foldername{i},'\','neuronIndividuals_new.mat']);
    for j=1:size(group,2)
        [~,regionn_debug{i,j},~,~,boundaryIdx{i,j},A_color_debug{i,j}]=cluster_region_cal_022422(group{i,j},neuronIndividuals_new,[]);
    end
end

% 2. find neighboring boundaries
all_nei_bd={};
all_nei_bd_CellIdx={};
all_nei_bd_closeIdx={};
all_nei_idx={};
for i=1:size(group,1)
    load([foldername{i},'\','neuronIndividuals_new.mat']);
    for j=1:size(group,2)
        bdCur=boundaryIdx{i,j};
        cen=neuron_centroid_calculation(neuronIndividuals_new{j},size(neuronIndividuals_new{j}.Cn));
        
        neighboor_dis_all=squareform(pdist(cen));
        neighboor_dis_all(neighboor_dis_all==0)=inf;
        neighboor_dis=quantile(min(neighboor_dis_all,[],2),0.95);

        % 1. extract all boundary fr eah cluster, and record their idx.
        all_bd={};
        all_bd_cellIdx={};
        all_bd_clustIdx={};
        ctt=1;
        for k1=1:length(bdCur)
            for k2=1:length(bdCur{k1})
                if ~isempty(bdCur{k1}{k2})
                    all_bd{ctt}=cen(bdCur{k1}{k2}(1:end-1),:);
                    all_bd_cellIdx{ctt}=bdCur{k1}{k2}(1:end-1);
                    all_bd_clustIdx{ctt}=ones(length(bdCur{k1}{k2}(1:end-1)),1)*k1;
                    ctt=ctt+1;
                end 
            end
        end
        
        % 2. neighboring clusters and index pairs  
        nei_bd={};
        nei_bd_CellIdx={};
        nei_bd_closeIdx={};
        nei_idx={};
        for k1=1:length(all_bd)-1 % each k1 and k2 is a individual cluster region for a specific cluster
            for k2=k1+1:length(all_bd)
                coor1=all_bd{k1};
                coor2=all_bd{k2};
                
                minDis=[];
                minDisIdx=[];
                for k3=1:size(coor1,1)
                    dis=sum((coor2-coor1(k3,:)).^2,2).^0.5;
                    minDis(k3)=min(dis);
                    if min(dis)<neighboor_dis
                        minDisIdx=[minDisIdx;[repmat(all_bd_cellIdx{k1}(k3),length(find(dis==min(dis))),1),all_bd_cellIdx{k2}(find(dis==min(dis)))]]; % just pick up the first one
                    end
                end
                % minDisIdx is the min distances of neurons in cluster
                % region k1 to neurons in cluster region k2
                if min(minDis)<neighboor_dis&&all_bd_clustIdx{k1}(1)~=all_bd_clustIdx{k2}(1) % at least one neuron pair
                    nei_idx=[nei_idx;{all_bd_clustIdx{k1},all_bd_clustIdx{k2}}];
                    nei_bd_CellIdx=[nei_bd_CellIdx;{all_bd_cellIdx{k1},all_bd_cellIdx{k2}}];
                    nei_bd_closeIdx=[nei_bd_closeIdx;{minDisIdx}];
                    nei_bd=[nei_bd;{coor1,coor2}];
                end
            end
        end
        
        all_nei_bd{i,j}=nei_bd;
        all_nei_bd_CellIdx{i,j}=nei_bd_CellIdx;
        all_nei_bd_closeIdx{i,j}=nei_bd_closeIdx;
        all_nei_idx{i,j}=nei_idx;
    end
end

% 3. cross correlation, intra-inter boundary neurons
all_crossCorr_intra={};
all_crossCorr_intra_border={};
all_crossCorr_inter_border={};
for i=1:size(group,1)
    load([foldername{i},'\','neuronIndividuals_new.mat']);
    for j=1:size(group,2)
        nC=neuronIndividuals_new{j}.C;
        nC=zscore(nC,[],2); % zscore to make the amplitude does not dominate the height of cross-corr
        cen=neuron_centroid_calculation(neuronIndividuals_new{j},size(neuronIndividuals_new{j}.Cn));

        neighboor_dis_all=squareform(pdist(cen));
        neighboor_dis_all(neighboor_dis_all==0)=inf;
        neighboor_dis=quantile(min(neighboor_dis_all,[],2),0.95);
        
        tic;
         % all cross correlation, all intra-cluster
        all_xcorr=[];
        ugp=unique(group{i,j});
        ugp(ugp==-1)=[];
        for k=1:length(ugp)
            nC1_all=nC(group{i,j}==ugp(k),:);
            cen_curGp=cen(group{i,j}==ugp(k),:);
            % attn: for better comparsion, we only pick intra-cluster pairs
            % within the neighbor_his
            
            ctt=1;
            for k1=1:size(nC1_all,1)-1
                for k2=k1+1:size(nC1_all,1)
                    cen1=cen_curGp(k1,:);
                    cen2=cen_curGp(k2,:);
                    if sum((cen1-cen2).^2)^0.5<=neighboor_dis % only pick the pairs with distance smaller than neighbor_dis
                        all_xcorr(ctt,:)=xcorr(nC1_all(k1,:),nC1_all(k2,:),maxLagg,'unbiased');
                        ctt=ctt+1;
                    end
                end
            end
        end
        
        % all cross correlation, borders
        all_xcorr_border=[];
        all_xcorr_border_close=[];
        for k=1:size(all_nei_bd_CellIdx{i,j},1)
   
            nC1_border=nC(unique(all_nei_bd_CellIdx{i,j}{k,1}),:);
            nC2_border=nC(unique(all_nei_bd_CellIdx{i,j}{k,2}),:);
            
            nC1_close=nC(unique(all_nei_bd_closeIdx{i,j}{k}(:,1)),:);
            nC2_close=nC(unique(all_nei_bd_closeIdx{i,j}{k}(:,2)),:);

            % all cross correlation, border intra-cluster
            ctt=1;
            for k1=1:size(nC1_border,1)-1
                for k2=k1+1:size(nC1_border,1)
                    cen_curGp1=cen(unique(all_nei_bd_CellIdx{i,j}{k,1}),:);
                    cen1=cen_curGp1(k1,:);
                    cen2=cen_curGp1(k2,:);
                    if sum((cen1-cen2).^2)^0.5<=neighboor_dis
                        all_xcorr_border(ctt,:)=xcorr(nC1_border(k1,:),nC1_border(k2,:),maxLagg,'unbiased');
                        ctt=ctt+1;
                    end
                end
            end
            
            for k1=1:size(nC2_border,1)-1
                for k2=k1+1:size(nC2_border,1)
                    cen_curGp2=cen(unique(all_nei_bd_CellIdx{i,j}{k,2}),:);
                    cen1=cen_curGp2(k1,:);
                    cen2=cen_curGp2(k2,:);
                    if sum((cen1-cen2).^2)^0.5<=neighboor_dis
                        all_xcorr_border(ctt,:)=xcorr(nC2_border(k1,:),nC2_border(k2,:),maxLagg,'unbiased');
                        ctt=ctt+1;
                    end
                end
            end
            
            % all cross correlation, border close inter-cluster
            ctt=1;
            for k1=1:size(nC1_close,1)
                for k2=1:size(nC2_close,1)
                    all_xcorr_border_close(ctt,:)=xcorr(nC1_close(k1,:),nC2_close(k2,:),maxLagg,'unbiased');
                    ctt=ctt+1;
                end
            end
            
        end
        toc;
        
        all_crossCorr_intra{i,j}=all_xcorr;
        all_crossCorr_intra_border{i,j}=all_xcorr_border;
        all_crossCorr_inter_border{i,j}=all_xcorr_border_close;
    end
end
