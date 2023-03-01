function [fr_all_laps_avg_sm]=smooth_laps_ratemaps(fr_all_laps_avg_line,h)

if isempty(h)
    h=[1]; % no conv
end

fr_all_laps_avg_sm={};
for tk=1:size(fr_all_laps_avg_line,1)
    for i=1:size(fr_all_laps_avg_line,2)       
        for j=1:size(fr_all_laps_avg_line{tk,i}{1},1)
            f1=fr_all_laps_avg_line{tk,i}{1}(j,:);
            f1=conv(f1,h,'same');
            fr_all_laps_avg_sm{tk,i}{1}(j,:)=f1;

            f2=fr_all_laps_avg_line{tk,i}{2}(j,:);
            f2=conv(f2,h,'same');
            fr_all_laps_avg_sm{tk,i}{2}(j,:)=f2;
        end
    end
end
