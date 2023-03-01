function trim_black_frames_finalConc(videoname,savname)

load(videoname);
sum1=nanmean(Y,1);
sum2=nanmean(squeeze(sum1),1);
idx=find(sum2<15);
if idx(1)==1
    idx=idx(2:end);
end
Y(:,:,idx)=[];
Ysiz=size(Y);
save(savname,'Y','Ysiz','-v7.3');