%suppl Fig11: length of transient, synchorny time
run('D:\Xu_clusterting_paper_prep11_2020\final_code\data_prepare\neuron_data_info.m')

%% 1. length of transient
% multiGeo
fname1={};
for i=1:length(foldername_multiGeo)
    fname1{i}=[foldername_multiGeo{i},'\','neuronIndividuals_new.mat'];
end

% Fig2
fname3={};
for i=1:length(foldername_fig2)
    fname3{i}=[foldername_fig2{i},'\','neuronIndividuals_new.mat'];
end

% multiGeo, transient length and (synchory time -- how)
[medianLength1,maxLength1,meanL1,idx1]=all_transient_duration(fname1,[3],15);
[medianLength3,maxLength3,meanL3,idx3]=all_transient_duration(fname3,[1],15);

%% MAX
subplot(4,1,1)
histogram(maxLength1,'binWidth',0.5)
line([quantile(maxLength1,0.99),quantile(maxLength1,0.99)],[0,800],'lineStyle','--','color','r')
xlim([0 55])

subplot(4,1,2)
histogram(maxLength3,'binWidth',0.5)
line([quantile(maxLength3,0.99),quantile(maxLength3,0.99)],[0,800],'lineStyle','--','color','r')
xlim([0 55])


set(gcf,'renderer','painters');

%% MEAN
subplot(4,1,1)
histogram(meanL1,'binWidth',0.25)
line([quantile(meanL1,0.99),quantile(meanL1,0.99)],[0,400],'lineStyle','--','color','r')
xlim([0 30])

subplot(4,1,2)
histogram(meanL3,'binWidth',0.25)
line([quantile(meanL3,0.99),quantile(meanL3,0.99)],[0,600],'lineStyle','--','color','r')
xlim([0 30])

set(gcf,'renderer','painters');

[N,EDGES]=histcounts(meanL1,'BinWidth',0.25)
[N,EDGES]=histcounts(meanL3,'BinWidth',0.25)

quantile(meanL1,0.99)
quantile(meanL3,0.99)


%% trace example
load(fname1{1})
plot(neuronIndividuals_new{4}.C(53,:));
hold on;

thres=3*std(neuronIndividuals_new{4}.S(53,:),[],2);
line([0,8900],[thres,thres],'lineStyle','--','color','r');

nCk=neuronIndividuals_new{4}.C(53,:)>thres;
stats=regionprops(nCk);

nCt=neuronIndividuals_new{4}.C(53,:);
for i=1:length(stats)
    if stats(i).Area/15>5
        % larger than 5sec area
        start=round(stats(i).BoundingBox(1));
        endd=round(stats(i).BoundingBox(1)+stats(i).BoundingBox(3));
        plot([start:endd],nCt(start:endd),'color','g');
    end
end

plot(neuronIndividuals_new{4}.S(53,:),'color','m');
