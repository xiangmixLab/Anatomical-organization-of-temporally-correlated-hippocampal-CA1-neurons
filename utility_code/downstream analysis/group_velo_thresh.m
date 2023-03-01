function [mean_thresh,median_thresh]=group_velo_thresh(all_velo)

 for i=1:size(all_velo,1)
    for j=1:size(all_velo,2)
        mean_velo(i,j)=nanmean(all_velo{i,j});
    end
 end

 mean_thresh=mean(mean_velo(:));
 median_thresh=median(mean_velo(:));
 