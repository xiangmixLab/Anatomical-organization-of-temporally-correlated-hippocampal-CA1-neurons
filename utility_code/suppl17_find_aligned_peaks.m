function [aligned_peaks,non_aligned_peaks,total_aligned_peaks,total_non_aligned_peaks,totalpeaks_nc1,totalpeals_nc2]=suppl17_find_aligned_peaks(nc1,nc2)

% suppose nc1 and nc2 are smoothed
aligned_peaks=cell(1,3);
non_aligned_peaks=cell(1,2);

total_aligned_peaks=0;
total_non_aligned_peaks=0;

[peak_locs1,transients1]=find_peaks_and_transients(nc1); % transient will be a array with [leftBound:1:rightBound]
[peak_locs2,transients2]=find_peaks_and_transients(nc2);

ctt=1;
for i=1:length(peak_locs1)
    leftSearchBound=max(peak_locs1(i)-3*15,1);
    rightSearchBound=min(peak_locs1(i)+3*15,length(nc2));
    
    allAlignedPeaksIdx=find((peak_locs2>=leftSearchBound).*(peak_locs2<=rightSearchBound)==1);
    
    % closest trim
    if length(allAlignedPeaksIdx)>1
         dis=abs(peak_locs2(allAlignedPeaksIdx)-peak_locs1(i));
         allAlignedPeaksIdx=allAlignedPeaksIdx(dis==min(dis));
    end
    
    % overlap trim
    if length(allAlignedPeaksIdx)>1
        overlaps=[];
        for j=1:length(allAlignedPeaksIdx)
            period1=transients1{i};
            period2=transients2{allAlignedPeaksIdx(j)};
            overlaps(j)=length(intersect(period1,period2));
        end
        
        allAlignedPeaksIdx=allAlignedPeaksIdx(find(overlaps==max(overlaps)));
        
        if length(allAlignedPeaksIdx)>1 % if there's still two or more peaks with the same overlap, bypass this peak.
            allAlignedPeaksIdx=[];
        else
            if max(overlaps)/length(period1)<0.5
                allAlignedPeaksIdx=[];
            end
        end
    end
    
    % if already picked, check distance, trim the previous record if the
    % previous distance is larger than current

    if ~isempty(aligned_peaks)&&~isempty(allAlignedPeaksIdx)
        if sum(cell2mat(aligned_peaks(:,2))==allAlignedPeaksIdx)>0
            idx=cell2mat(aligned_peaks(:,2))==allAlignedPeaksIdx;
            dis_prev=aligned_peaks{idx,3};
            dis_curr=abs(peak_locs1(i)-peak_locs2(allAlignedPeaksIdx));
            if dis_curr<dis_prev
                aligned_peaks(idx,:)=[];
            end
        end
    end
    
    if length(allAlignedPeaksIdx)>0
        aligned_peaks{ctt,1}=peak_locs1(i);
        aligned_peaks{ctt,2}=peak_locs2(allAlignedPeaksIdx);
        aligned_peaks{ctt,3}=abs(peak_locs1(i)-peak_locs2(allAlignedPeaksIdx)); % distance
        ctt=ctt+1;
    end
end

aligned_peaks1=cell2mat(aligned_peaks(:,1));
aligned_peaks2=cell2mat(aligned_peaks(:,2));

non_align_peaks1=setdiff(peak_locs1,aligned_peaks1);
non_align_peaks2=setdiff(peak_locs2,aligned_peaks2);

non_aligned_peaks{1}=non_align_peaks1;
non_aligned_peaks{2}=non_align_peaks2;

total_aligned_peaks=length(aligned_peaks1);
total_non_aligned_peaks=length(non_align_peaks1)+length(non_align_peaks2);

totalpeaks_nc1=length(peak_locs1);
totalpeals_nc2=length(peak_locs2);

% check
% plot(nc1);hold on;plot(nc2);
% plot(aligned_peaks1,ones(length(aligned_peaks1),1)*5,'*');
% plot(non_align_peaks1,ones(length(non_align_peaks1),1)*5,'^');
% plot(non_align_peaks2,ones(length(non_align_peaks2),1)*5,'^');