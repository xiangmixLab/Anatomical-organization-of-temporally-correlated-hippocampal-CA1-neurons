function rm_line=linearTrack_ratemap2line(fr)

rm_line={};
for i=1:size(fr,1)
    for j=1:size(fr,2)
        frt=[];
        
        maxsize=[];
        for k=1:length(fr{i,j})
            if ~isempty(fr{i,j}{k})
                maxsize=[maxsize;size(fr{i,j}{k})];
            end
        end
        maxsize=max(maxsize,[],1);
        
        for k=1:length(fr{i,j}) 
            if ~isempty(fr{i,j}{k})
                if size(fr{i,j}{k},1)<size(fr{i,j}{k},2)
                    frt(k,:)=nanmean(fr{i,j}{k},1);
                else
                    frt(k,:)=nanmean(fr{i,j}{k},2);
                end
            else
                if maxsize(1)<maxsize(2) % horizontal
                    frt(k,:)=zeros(1,maxsize(2));
                else % vertical
                    frt(k,:)=zeros(1,maxsize(1));
                end
            end
        end
        rm_line{i,j}=frt;
    end
end