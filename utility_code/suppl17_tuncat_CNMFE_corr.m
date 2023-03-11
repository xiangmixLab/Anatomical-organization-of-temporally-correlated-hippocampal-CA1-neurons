function [corr_all,corr_per_mice]=suppl17_tuncat_CNMFE_corr(dpath)

figure;
corr_all=[];
corr_per_mice={};
ctt=1;
for i=1:length(dpath)
    load([dpath{i},'\','neuronIndividuals_new.mat'])
    load([dpath{i},'\','neuronIndividuals_new_tuncat.mat'])
    for j=1:length(neuronIndividuals_new)
        corr_curr=[];
        nC1=zscore(neuronIndividuals_new{1}.C,[],2);
        nC2=zscore(neuronIndividuals_new_tuncat{1}.C,[],2); % normalize amplitude difference
        for k=1:size(neuronIndividuals_new{j}.C,1)
            
            nc1=nC1(k,:)-min(nC1(k,:));
            nc2=nC2(k,:)-min(nC2(k,:));
            
            nc1=smoothdata(nc1,'gaussian',6*15);
            nc2=smoothdata(nc2,'gaussian',6*15);
            
            Cpeaks1=C_to_peakS(nc1);
            Cpeaks2=C_to_peakS(nc2);
            
            thresh1=0.1*max(Cpeaks1);
            thresh2=0.1*max(Cpeaks2);
            
            nc1(nc1<thresh1)=0;
            nc2(nc2<thresh2)=0;
            
            c=corrcoef(nc1,nc2);
            corr_all(ctt,1)=c(2);
            corr_curr(k,1)=c(2);
            ctt=ctt+1;
        end
        corr_per_mice{i,j}=corr_curr;
    end
end
histogram(corr_all);
set(gcf,'renderer','painters');
