function [locs,transients]=find_peaks_and_transients(nc1)

[~,locs]=findpeaks(nc1);
transients={};

for i=1:length(locs)
    curr_locs=locs(i);
    % left search
    leftBound=curr_locs;
    while nc1(leftBound)>0 && leftBound>1
        leftBound=leftBound-1;
    end
    
    % right search
    rightBound=curr_locs;
    while nc1(rightBound)>0 && rightBound<length(nc1)
        rightBound=rightBound+1;
    end
    transients{i}=[leftBound:rightBound];
end
