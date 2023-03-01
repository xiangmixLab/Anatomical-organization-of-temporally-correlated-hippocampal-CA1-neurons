function [fr_maxf,fr_bin_maxf]=max_rate_field_select(fr)

for i=1:length(fr)
    frt=fr{i}>0.4*max(fr{i}(:));
    
    stats=regionprops(frt);
    
    peak_rates=[];
    for j=1:length(stats)
        bbox=floor(stats(j).BoundingBox);
        bbox(bbox==0)=1;
        rate_seg=frt(bbox(2):bbox(2)+bbox(4)-1,bbox(1):bbox(1)+bbox(3)-1);
        peak_rates=max(rate_seg(:));
    end
    
    field_select=find(peak_rates==max(peak_rates));
    
    fr_maxf{i,1}=frt*0;
    bbox=floor(stats(field_select).BoundingBox);
    bbox(bbox==0)=1;
    fr_maxf{i,1}(bbox(2):bbox(2)+bbox(4)-1,bbox(1):bbox(1)+bbox(3)-1)=frt(bbox(2):bbox(2)+bbox(4)-1,bbox(1):bbox(1)+bbox(3)-1);
    
    fr_bin_maxf{i,1}=fr_maxf{i,1}>0;
end
