
function [all_dis_corr_circle_all_intra,all_dis_corr_circle_all_inter,midx,f1,f2,gof1,gof2]=pairwise_dis_tempCorr_071421(neuron_all,group_all)
all_dis_circle_intra={};
all_corr_circle_intra={};
all_dis_circle_inter={};
all_corr_circle_inter={};

midx={};
midx_intra=[];
midx_inter=[];

for i=1:size(neuron_all,1)
    for j=1:size(neuron_all,2)
        group=group_all{i,j};
        uni_group=unique(group);

        dis_all=squareform(pdist(neuron_all{i,j}.centroid))*2; % *2 is added 042320, in light of the resampling of the frame make all centroid distance only half of what it should be
        corr_all=squareform(1-pdist(neuron_all{i,j}.C,'correlation'));

        % in group
        tic;
        for k1=1:length(uni_group)
            %distance
            t=triu(dis_all(group==k1,group==k1),1);
            mask = tril(true(size(corr_all(group==k1,group==k1))),1);
            all_dis_circle_intra{i,j}{k1}=t(~mask);
            t=triu(dis_all(group==k1,group~=k1),1);
            mask = tril(true(size(corr_all(group==k1,group~=k1))),1);
            all_dis_circle_inter{i,j}{k1}=t(~mask);
            %correlation
            t=triu(corr_all(group==k1,group==k1),1);
            mask = tril(true(size(corr_all(group==k1,group==k1))),1);
            all_corr_circle_intra{i,j}{k1}=t(~mask);
            t=triu(corr_all(group==k1,group~=k1),1);
            mask = tril(true(size(corr_all(group==k1,group~=k1))),1);
            all_corr_circle_inter{i,j}{k1}=t(~mask);   
            toc;
        end
    end
end

all_dis_corr_circle_all_intra=[];
all_dis_corr_circle_all_inter=[];
ctt1=1;
ctt2=1;
for k=1:size(neuron_all,1)
    for j=1:size(neuron_all,2)  
        group=group_all{k,j};
        uni_group=unique(group);

        for i=1:length(unique(group))
            all_dis_corr_circle_all_intra(ctt1:ctt1+length(all_dis_circle_intra{k,j}{i})-1,1)=all_dis_circle_intra{k,j}{i};
            all_dis_corr_circle_all_intra(ctt1:ctt1+length(all_dis_circle_intra{k,j}{i})-1,2)=all_corr_circle_intra{k,j}{i};
            midx_intra(ctt1:ctt1+length(all_dis_circle_intra{k,j}{i})-1)=k;
            ctt1=ctt1+length(all_dis_circle_intra{k,j}{i});

            all_dis_corr_circle_all_inter(ctt2:ctt2+length(all_dis_circle_inter{k,j}{i})-1,1)=all_dis_circle_inter{k,j}{i};
            all_dis_corr_circle_all_inter(ctt2:ctt2+length(all_dis_circle_inter{k,j}{i})-1,2)=all_corr_circle_inter{k,j}{i};
            midx_inter(ctt2:ctt2+length(all_dis_circle_inter{k,j}{i})-1)=k;
            ctt2=ctt2+length(all_dis_circle_inter{k,j}{i});
        end
    end
end

midx={midx_intra,midx_inter};

non_nan_idx_intra=logical(~isnan(all_dis_corr_circle_all_intra(1:end,1)).*~isnan(all_dis_corr_circle_all_intra(1:end,2)));
non_nan_idx_inter=logical(~isnan(all_dis_corr_circle_all_inter(1:end,1)).*~isnan(all_dis_corr_circle_all_inter(1:end,2)));

[f1,gof1]=fit(all_dis_corr_circle_all_intra(non_nan_idx_intra,1),all_dis_corr_circle_all_intra(non_nan_idx_intra,2),'power2');
[f2,gof2]=fit(all_dis_corr_circle_all_inter(non_nan_idx_inter,1),all_dis_corr_circle_all_inter(non_nan_idx_inter,2),'power2');

figure;
hold on;
h1=plot(f1,all_dis_corr_circle_all_intra(1:40:end,1),all_dis_corr_circle_all_intra(1:40:end,2));
h2=plot(f2,all_dis_corr_circle_all_inter(1:40:end,1),all_dis_corr_circle_all_inter(1:40:end,2));

set([h1],'color',[0.9290    0.6940    0.1250],'MarkerSize',10)
set([h2],'color',[0.4940    0.1840    0.5560],'MarkerSize',10)
