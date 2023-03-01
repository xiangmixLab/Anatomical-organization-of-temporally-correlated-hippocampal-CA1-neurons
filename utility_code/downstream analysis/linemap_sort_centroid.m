function [idxx]=linemap_sort_centroid(fr_bin_11,fr_sm_11)

for i=1:length(fr_bin_11)
    fr1=fr_bin_11{i}; 
    fs1=fr_sm_11{i};
    
    cen=[];
    for j=1:size(fr1,1)
        stats=regionprops(fr1(j,:));
        if ~isempty(stats)
            cent=[];
            intensityt=[];
            for k=1:length(stats)
                intensityt=[intensityt;sum(fs1(j,round(stats(k).BoundingBox(2))+1:round(stats(k).BoundingBox(2)+stats(k).BoundingBox(4))))];
                cent=[cent;stats(k).Centroid];
            end
            cen=[cen;cent(intensityt==max(intensityt),:)];
        else
            cen=[cen;[-1 -1]];
        end
    end
    
    if max(abs(cen(:,1)))<=2
        [~,idxx{i}]=sort(cen(:,2));
    else
        [~,idxx{i}]=sort(cen(:,1));     
    end
end
    