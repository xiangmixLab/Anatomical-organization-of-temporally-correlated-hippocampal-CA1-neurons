function [fr_all_laps_avg_sm]=smooth_linemap_ratemaps(fr_all_laps_avg_line,h)

if isempty(h)
    h=[1]; % no conv
end

fr_all_laps_avg_sm={};
for tk=1:size(fr_all_laps_avg_line,1)
    for i=1:size(fr_all_laps_avg_line,2)       
        for j=1:size(fr_all_laps_avg_line{tk,i},1)
            f1=fr_all_laps_avg_line{tk,i}(j,:);
            f1=conv(f1,h,'same');
            fr_all_laps_avg_sm{tk,i}(j,:)=f1;
        end
    end
end
