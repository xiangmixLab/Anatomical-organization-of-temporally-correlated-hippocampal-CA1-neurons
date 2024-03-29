%% supplemental Fig.7: segregation analysis
run('D:\Xu_clusterting_paper_prep11_2020\final_code\data_prepare\neuron_data_info.m')
all_colors=distinguishable_colors(10);
load(['D:\Xu_clusterting_paper_prep11_2020\final_code\final_cluster_data\cluster_optimal_num\Fig2_cir_rec_clust_original.mat'])

%% data prepare for following analysis
%cluster
group_dat=group_ori_fig2;
foldername=foldername_fig2;

% behav
all_behav={};
for i=1:length(foldername)
    load([foldername{i},'\','Behav.mat']);
    all_behav{i,1}=behavIndividuals{1};
    all_behav{i,2}=behavIndividuals{2};
end

% ensemble firing
nnew_all={};
n_gp={};
n_gp_ensem={};
overlap_idx1={};
for tk=1:12
    load([foldername{tk},'\','neuronIndividuals_new.mat'])
    for l=1:2
        for i=1:length(unique(group_dat{tk,l}))
            n_gp{tk,l}{i}=current_gp_neuron(neuronIndividuals_new{1},group_dat{tk,l},i);
            
            n_gp1_C=n_gp{tk,l}{i}.C;
            n_gp1_S=C_to_peakS(n_gp1_C);
            n_gp1_C_ensem=nanmean(zscore(n_gp1_C,[],2),1);
            n_gp1_C_ensem(n_gp1_C_ensem<0)=0;
            n_gp1_S_ensem=C_to_peakS(n_gp1_C_ensem);
    
            n_gp_ensem{tk,l}.A=n_gp{tk,l}{i}.A;
            n_gp_ensem{tk,l}.C(i,:)=n_gp1_C_ensem;
            n_gp_ensem{tk,l}.C_raw(i,:)=n_gp1_C_ensem;
            n_gp_ensem{tk,l}.S(i,:)=n_gp1_S_ensem;
            n_gp_ensem{tk,l}.time=n_gp{tk,l}{i}.time;
            
        end
    end
end

% ensemble rate map
binsize=10;
firingrateAll_ensem_all={};
countTime_ensem_all={};
for i=1:size(group_dat,1)
    for j=1:size(group_dat,2)
        uni_group=unique(group_dat{i,j});
        pct=1;
        [firingrateAll_ensem,~,~,countTime_ensem] = calculatingCellSpatialForSingleData_040321(n_gp_ensem{i,j},all_behav{i,j}.position,all_behav{i,j}.time,[0,0,max(all_behav{i,j}.position,[],1)],binsize,1:size(n_gp_ensem{i,j}.C,1),3*std(n_gp_ensem{i,j}.C,[],2),'S',[],[],[0 inf],10);
     
        firingrateAll_ensem_all{i,j}=firingrateAll_ensem;
        countTime_ensem_all{i,j}=countTime_ensem;
        
    end
end

%% overlap calculation
overlap_idx1={};
for tk=1:size(group_dat,1)
    for k=1:size(group_dat,2)
        overlap_idx=[];
        ctt=1;
        for t=1:size(firingrateAll_ensem_all{tk,k},2)-1
            for j=t+1:size(firingrateAll_ensem_all{tk,k},2)
                rt1=firingrateAll_ensem_all{tk,k}{t};
                rt2=firingrateAll_ensem_all{tk,k}{j};
                rt1=filter2DMatrices(rt1,1);
                rt2=filter2DMatrices(rt2,1);
                rt1_bin=rt1>max(rt1(:))*0.5;
                rt2_bin=rt2>max(rt2(:))*0.5;
                overlap_idx(ctt)=sum(sum(rt1_bin.*rt2_bin))*2/(sum(sum((rt1_bin+rt2_bin)>0)));
                ctt=ctt+1;
            end
        end
        
        overlap_idx1{tk,k}=overlap_idx;
    end
end

%% ensemble ratemap overlap
% 50% field example
fr1=firingrateAll_ensem_all{2,1}{1};
fr2=firingrateAll_ensem_all{9,2}{1};
ctime1=countTime_ensem_all{2,1};
ctime2=countTime_ensem_all{9,2};

[bx1,by1]=boundary_from_binImg(fr1,0.5,1);
subplot(2,1,1);ratemap_plot(fr1,ctime1,1,0,[]);hold on; plot(bx1{1},by1{1},'--','color','k'); % Ntg old
[bx2,by2]=boundary_from_binImg(fr2,0.5,1);
subplot(2,1,2);ratemap_plot(fr2,ctime2,1,0,[]);hold on; plot(bx2{1},by2{1},'--','color','k'); % Ntg old

colormap(jet)

% illustrate overlap
binsize=10;

figure;
colorclust=distinguishable_colors(10);
colorclust(4,:)=colorclust(10,:)
rcomb_all={};
for tk=1:size(group_dat,1)
    for k=1:size(group_dat,2)
        rcomb=zeros(length(0:binsize:ceil(max(all_behav{tk,k}.position(:,2)))),length(0:binsize:ceil(max(all_behav{tk,k}.position(:,1)))),3);
        uni_gp=unique(group_dat{tk,k});

        
        for i=1:length(uni_gp)
% 
            rs=filter2DMatrices(firingrateAll_ensem_all{tk,k}{i},1);
            r1=rs>0.5*max(rs(:));
            [ridxi,ridxj]=find(r1==1);

            for k1=1:length(ridxi)
                if sum(rcomb(ridxi(k1),ridxj(k1),:))==0
                    rcomb(ridxi(k1),ridxj(k1),:)=colorclust(i,:);
                else
                    rcomb(ridxi(k1),ridxj(k1),:)=(squeeze(rcomb(ridxi(k1),ridxj(k1),:))'+colorclust(i,:))./2;
                end
                
            end
            
            for k1=1:size(rcomb,3)
                rt=squeeze(rcomb(:,:,k1));
                rt(countTime_ensem_all{tk,k}==0)=255;
                rcomb(:,:,k1)=rt;
            end
        end
        rcomb_all{tk,k}=rcomb;
        subplot(size(group_dat,2),size(group_dat,1),tk+(k-1)*size(group_dat,1))
        imagesc(rcomb_all{tk,k});
        axis ij
        shading flat;
        axis image
        axis off
    end
end

colorbk=zeros(1,10,3);
for i=1:size(colorclust,1)
    colorbk(1,i,:)=colorclust(i,:);
end
imagesc(colorbk)
axis off

% illustrate overlap distribution
all_overlap_cir=[];
all_overlap_sqr=[];
for i=1:12
    
    all_overlap_cir=[all_overlap_cir,overlap_idx1{i,1}];
    all_overlap_sqr=[all_overlap_sqr,overlap_idx1{i,2}];
end
subplot(121);
histogram(all_overlap_cir,'binWidth',0.1); hold on;
ptile79=quantile(all_overlap_cir,0.86);
plot(ones(1,100)*ptile79,[1:100],'--')
set(gcf,'renderer','painters');
subplot(122);
histogram(all_overlap_sqr,'binWidth',0.1); hold on;
ptile67=quantile(all_overlap_sqr,0.75);
plot(ones(1,100)*ptile67,[1:100],'--')
set(gcf,'renderer','painters');

