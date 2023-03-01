function [rm_line_dir1,rm_line_dir2]=linearTrack_ratemap2line_laps(fr)

rm_line={};
for i=1:size(fr,1)
    for j=1:size(fr,2)
        fr1=fr{i,j};
        for k1=1:size(fr1,1)
            for k2=1:size(fr1,2)
                if ~isempty(fr1)
                frt=[];

                maxsize=[];
                for k=1:length(fr1{k1,k2})
                    if ~isempty(fr1{k1,k2}{k})
                        maxsize=[maxsize;size(fr1{k1,k2}{k})];
                    end
                end
                maxsize=max(maxsize,[],1);

                for k=1:length(fr1{k1,k2}) 
                    if ~isempty(fr1{k1,k2}{k})
                        if size(fr1{k1,k2}{k},1)<size(fr1{k1,k2}{k},2)
                            frt(k,:)=nanmean(fr1{k1,k2}{k},1);
                        else
                            frt(k,:)=nanmean(fr1{k1,k2}{k},2);
                        end
                    else
                        if maxsize(1)<maxsize(2) % horizontal
                            frt(k,:)=zeros(1,maxsize(2));
                        else % vertical
                            frt(k,:)=zeros(1,maxsize(1));
                        end
                    end
                end
                rm_line{i,j}{k1,k2}=frt;
                end
            end
        end
    end
end

