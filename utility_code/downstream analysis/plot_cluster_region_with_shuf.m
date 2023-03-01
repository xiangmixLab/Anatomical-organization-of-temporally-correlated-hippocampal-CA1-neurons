function plot_cluster_region_with_shuf(avg_region,avg_region_shuf,yrange)

avg_region_cond1_mat=[];
for i=1:size(avg_region,1)
    avg_region_cond1_mat(i,:)=cell2mat(avg_region{i,1}); 
end

avg_region_cond1_shuf_mat=[];
for i=1:size(avg_region_shuf,1) 
    shuf_region=avg_region_shuf{i,1};
    shuf_region_mat=[];
    for j=1:length(shuf_region)
        shuf_region_mat(j)=mean(shuf_region{j});
    end
    avg_region_cond1_shuf_mat(i,:)=shuf_region_mat;
end

avg_region_cond1_shuf_95_mat=[];
for i=1:size(avg_region_shuf,1) 
    shuf_region=avg_region_shuf{i,1};
    shuf_region_mat=[];
    for j=1:length(shuf_region)
        shuf_region_mat(j)=quantile(shuf_region{j},0.95);
    end
    avg_region_cond1_shuf_95_mat(i,:)=shuf_region_mat;
end
% 
% for i=1:5
%     subplot(1,5,i)
%     plot(avg_region_cond1_mat(i,:));hold on;
%     plot(avg_region_cond1_shuf_95_mat(i,2:end));
%     plot(avg_region_cond1_shuf_mat(i,2:end));
% end

% 840F neuron density is smaller than the other four
errorbar(mean(avg_region_cond1_mat,1),std(avg_region_cond1_mat,[],1)/(size(avg_region_cond1_mat,1)^0.5),'-o','color',[0    0.4470    0.7410])
hold on;
% errorbar(mean(avg_region_cond1_shuf_95_mat(:,2:end),1),std(avg_region_cond1_shuf_95_mat(:,2:end),[],1)/(size(avg_region_cond1_shuf_95_mat(:,2:end),1)^0.5),'-s','color',[0.8500    0.3250    0.0980])
errorbar(mean(avg_region_cond1_shuf_mat(:,2:end),1),std(avg_region_cond1_shuf_mat(:,2:end),[],1)/(size(avg_region_cond1_shuf_mat(:,2:end),1)^0.5),'-^','color',[0.9290    0.6940    0.1250])
ylim(yrange);