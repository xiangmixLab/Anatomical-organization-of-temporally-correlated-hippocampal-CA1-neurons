function [peakDiff,peakDiff_percentage]=suppl17_tuncat_CNMFE_peakDiff(dpath)

%%%%%%%%%% # of peaks above threshold difference 
peakDiff=[];
peakDiff_percentage=[];
ctt=1;
for i=1:6
    load([dpath{i},'\','neuronIndividuals_new_cleaned.mat'])
    load([dpath{i},'\','neuronIndividuals_new_tuncat_cleaned.mat'])
    for j=1:1
        nC1=zscore(neuronIndividuals_new{1}.C,[],2);
        nC2=zscore(neuronIndividuals_new_tuncat{1}.C,[],2); % normalize amplitude difference
        
        for k=1:size(neuronIndividuals_new{j}.C,1)
            
            nc1=nC1(k,:)-min(nC1(k,:));
            nc2=nC2(k,:)-min(nC2(k,:));
            
            nc1=smoothdata(nc1,'gaussian',6*15);
            nc2=smoothdata(nc2,'gaussian',6*15);
            
            Cpeaks1=C_to_peakS(nc1);
            Cpeaks2=C_to_peakS(nc2);
            
            thresh1=0.2*max(Cpeaks1);
            thresh2=0.2*max(Cpeaks2);
            
            nc1(nc1<thresh1)=0;
            nc2(nc2<thresh2)=0;
            
            [aligned_peaks,non_aligned_peaks,total_aligned_peaks,total_non_aligned_peaks,nc1_pk,nc2_pk]=suppl17_find_aligned_peaks(nc1,nc2);

            peakDiff(ctt,1)=total_non_aligned_peaks;
            peakDiff_percentage(ctt,1)=total_non_aligned_peaks/(nc1_pk+nc2_pk);

            plot(nC1(k,:));hold on;plot(nC2(k,:)-10);
            line([0,size(nC1(k,:),2)],[thresh1,thresh1],'lineStyle','--','color','r');
            line([0,size(nC2(k,:),2)],[thresh2,thresh2]-10,'lineStyle','--','color','r');

            plot(cell2mat(aligned_peaks(:,1)),ones(length(cell2mat(aligned_peaks(:,1))),1)*5,'*','color','g');
            plot(non_aligned_peaks{1},ones(length(non_aligned_peaks{1}),1)*20,'^','color','r');
            plot(non_aligned_peaks{2},ones(length(non_aligned_peaks{2}),1)*20,'^','color','r');

            set(gcf,'renderer','painters');
            saveas(gcf,['D:\Xu_clusterting_paper_prep11_2020\Round26\suppl17_panels\peak_alignment\',num2str(i),'_',num2str(j),'_',num2str(k),'.eps'],'epsc');
            saveas(gcf,['D:\Xu_clusterting_paper_prep11_2020\Round26\suppl17_panels\peak_alignment\',num2str(i),'_',num2str(j),'_',num2str(k),'.png']);
            ctt=ctt+1;
            close;

        end
    end
end
figure;
histogram(peakDiff);
set(gcf,'renderer','painters');
figure;
histogram(peakDiff_percentage);
set(gcf,'renderer','painters');
