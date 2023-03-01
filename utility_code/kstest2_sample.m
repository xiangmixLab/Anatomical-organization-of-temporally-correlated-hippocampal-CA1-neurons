function mean_p=kstest2_sample(dat1,dat2,sampling_percentile)

p=[];

parfor i=1:1000
    idx1=randperm(length(dat1));
    sampleIdx1=idx1(1:floor(length(dat1)*sampling_percentile));
    idx2=randperm(length(dat2));
    sampleIdx2=idx2(1:floor(length(dat2)*sampling_percentile));
    
    dat1_part=dat1(sampleIdx1);
    dat2_part=dat2(sampleIdx2);
    
    [~,p(i)]=kstest2(dat1_part,dat2_part);
end

mean_p=mean(p);
    



