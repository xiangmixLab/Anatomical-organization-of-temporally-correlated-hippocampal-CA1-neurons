function cluster_crossCorr_illustration(group,foldername,miceIdx,neiIdx)

% 1. calculating boundary idx
for i=miceIdx
    load([foldername{i},'\','neuronIndividuals_new.mat']);
    for j=1:1
        [~,regionn_debug{i,j},~,~,boundaryIdx{i,j},A_color_debug{i,j}]=cluster_region_cal_022422(group{i,j},neuronIndividuals_new,[]);
    end
end

% 2. find neighboring boundaries
all_nei_bd={};
all_nei_bd_CellIdx={};
all_nei_bd_closeIdx={};
all_nei_idx={};
for i=miceIdx
    load([foldername{i},'\','neuronIndividuals_new.mat']);
    for j=1:1
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
                    if min(dis)<=neighboor_dis
                        minDisIdx=[minDisIdx;[repmat(all_bd_cellIdx{k1}(k3),length(find(dis==min(dis))),1),all_bd_cellIdx{k2}(find(dis==min(dis)))]]; % just pick up the first one
                    end
                end
                % minDisIdx is the min distances of neurons in cluster
                % region k1 to neurons in cluster region k2
                if min(minDis)<=neighboor_dis&&all_bd_clustIdx{k1}(1)~=all_bd_clustIdx{k2}(1) % at least one neuron pair
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

% test bound idx
% imagesc(A_color_debug{miceIdx,1})
% for i=1:length(all_bd)
%      text(all_bd{i}(:,1),all_bd{i}(:,2),num2str(i));hold on;
% end
% 
% imagesc(A_color_debug{miceIdx,1})
% text(all_nei_bd{miceIdx,1}{2,1}(:,1),all_nei_bd{miceIdx,1}{2,1}(:,2),num2str(i));hold on;

% 3. illustration
load([foldername{miceIdx},'\','neuronIndividuals_new.mat']);
cen=neuron_centroid_calculation(neuronIndividuals_new{1},size(neuronIndividuals_new{1}.Cn));

neighboor_dis_all=squareform(pdist(cen));
neighboor_dis_all(neighboor_dis_all==0)=inf;
neighboor_dis=quantile(min(neighboor_dis_all,[],2),0.95);

subplot(131) % all intra
imagesc(A_color_debug{miceIdx,1});
hold on;

idx=find(group{miceIdx,1}==1);
dist=squareform(pdist(cen(idx,:)));
dist(dist==0)=inf;
idx(min(dist,[],2)>neighboor_dis)=[];
plot(cen(idx,1),cen(idx,2),'*','markerSize',10,'color','y');

idx=find(group{miceIdx,1}==2);
dist=squareform(pdist(cen(idx,:)));
dist(dist==0)=inf;
idx(min(dist,[],2)>neighboor_dis)=[];
plot(cen(idx,1),cen(idx,2),'*','markerSize',10,'color','k');
% 
subplot(132) % all intra
imagesc(A_color_debug{miceIdx,1});
hold on;

idx=all_nei_bd_CellIdx{miceIdx,1}{neiIdx,1};
dist=squareform(pdist(cen(idx,:)));
dist(dist==0)=inf;
idx(min(dist,[],2)>neighboor_dis)=[];
plot(cen(idx,1),cen(idx,2),'*','markerSize',10,'color','y');

idx=all_nei_bd_CellIdx{miceIdx,1}{neiIdx,2};
dist=squareform(pdist(cen(idx,:)));
dist(dist==0)=inf;
idx(min(dist,[],2)>neighboor_dis)=[];
plot(cen(idx,1),cen(idx,2),'*','markerSize',10,'color','k');

subplot(133) % all intra
imagesc(A_color_debug{miceIdx,1});
hold on;
plot(cen(unique(all_nei_bd_closeIdx{miceIdx,1}{neiIdx}(:,1)),1),cen(unique(all_nei_bd_closeIdx{miceIdx,1}{neiIdx}(:,1)),2),'*','markerSize',10,'color','y');
plot(cen(unique(all_nei_bd_closeIdx{miceIdx,1}{neiIdx}(:,2)),1),cen(unique(all_nei_bd_closeIdx{miceIdx,1}{neiIdx}(:,2)),2),'*','markerSize',10,'color','k');
set(gcf,'renderer','painters')