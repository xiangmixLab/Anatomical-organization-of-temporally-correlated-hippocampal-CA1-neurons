function trim_black_frames(videoname)

load(videoname);
sum1=sum(concatenatedvideo_res,1);
sum2=sum(squeeze(sum1),1);
idx=find(sum2==0);
if idx(1)==1
    idx=idx(2:end);
end
concatenatedvideo_res(:,:,idx)=[];
save(videoname,'concatenatedvideo_res','-v7.3');