function [idxx]=laps_sort_centroid(fr_bin_11)

for i=1:length(fr_bin_11)
    fr1=fr_bin_11{i}{1}; % right
    fr2=fr_bin_11{i}{2}; % left
    
    cen={[],[]};
    for j=1:size(fr1,1)
        stats=regionprops(fr1(j,:));
        if ~isempty(stats)
            cen{1}=[cen{1};stats.Centroid];
        else
            cen{1}=[cen{1};[-1 -1]];
        end
        stats=regionprops(fr2(j,:));
        if ~isempty(stats)
            cen{2}=[cen{2};stats.Centroid];
        else
            cen{2}=[cen{2};[-1 -1]];
        end
    end
    
    if sum(diff(cen{1}(:,1)))==0
        [~,idxx{i}{1}]=sort(cen{1}(:,2));
        [~,idxx{i}{2}]=sort(cen{2}(:,2));
    else
        [~,idxx{i}{1}]=sort(cen{1}(:,1));
        [~,idxx{i}{2}]=sort(cen{2}(:,1));      
    end
end
    