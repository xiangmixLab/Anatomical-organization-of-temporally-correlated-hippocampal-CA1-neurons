function [fr_all_dir,ct_all_dir,ctime_all_dir]=laps_avg_fr_ct(fr,ct,ctime,nnum)
%% fr
fr_all_dir={};
for i=1:size(fr,1)
    for j=1:size(fr,2)
        fr1=[];
        fr2=[];
        ctt1=ones(nnum(i),1);
        ctt2=ones(nnum(i),1);
        for k1=1:size(fr{i,j},1)
            if ~isempty(fr{i,j}{k1,1})
                if isempty(fr1)
                    fr1=fr{i,j}{k1,1};
                else
                    for k=1:length(fr{i,j}{k1,1})
                        if isempty(fr1{k})&&~isempty(fr{i,j}{k1,1}{k})
                            fr1{k}=fr{i,j}{k1,1}{k};
                        end
                        if ~isempty(fr1{k})&&isempty(fr{i,j}{k1,1}{k})
                            fr1{k}=fr1{k}+fr1{k}*0;
                        end  
                        if ~isempty(fr1{k})&&~isempty(fr{i,j}{k1,1}{k})
                            fr1{k}=fr1{k}+fr{i,j}{k1,1}{k};
                            ctt1(k)=ctt1(k)+1;
                        end 
                    end
                end
                
            end
            if ~isempty(fr{i,j}{k1,2})
                if isempty(fr2)
                    fr2=fr{i,j}{k1,2};
                else
                    for k=1:length(fr{i,j}{k1,2})
                        if isempty(fr2{k})&&~isempty(fr{i,j}{k1,2}{k})
                            fr2{k}=fr{i,j}{k1,2}{k};
                        end
                        if ~isempty(fr2{k})&&isempty(fr{i,j}{k1,2}{k})
                            fr2{k}=fr2{k}+fr2{k}*0;
                        end  
                        if ~isempty(fr2{k})&&~isempty(fr{i,j}{k1,2}{k})
                            fr2{k}=fr2{k}+fr{i,j}{k1,2}{k};
                            ctt2(k)=ctt2(k)+1;
                        end 
                    end
                end
            end
        end
        
        for k=1:length(fr1)
            fr1{k}=fr1{k}/ctt1(k);
        end
        for k=1:length(fr2)
            fr2{k}=fr2{k}/ctt2(k);
        end
        
        fr_all_dir{1}{i,j}=fr1;
        fr_all_dir{2}{i,j}=fr2;
    end
end
 %% ct
ct_all_dir={};
for i=1:size(ct,1)
    for j=1:size(ct,2)
        ct1=[];
        ct2=[];
        ctt1=ones(nnum(i),1);
        ctt2=ones(nnum(i),1);
        for k1=1:size(ct{i,j},1)
            if ~isempty(ct{i,j}{k1,1})
                if isempty(ct1)
                    ct1=ct{i,j}{k1,1};
                else
                    for k=1:length(ct{i,j}{k1,1})
                        if isempty(ct1{k})&&~isempty(ct{i,j}{k1,1}{k})
                            ct1{k}=ct{i,j}{k1,1}{k};
                        end
                        if ~isempty(ct1{k})&&isempty(ct{i,j}{k1,1}{k})
                            ct1{k}=ct1{k}+ct1{k}*0;
                        end  
                        if ~isempty(ct1{k})&&~isempty(ct{i,j}{k1,1}{k})
                            ct1{k}=ct1{k}+ct{i,j}{k1,1}{k};
                            ctt1(k)=ctt1(k)+1;
                        end 
                    end
                end
                
            end
            if ~isempty(ct{i,j}{k1,2})
                if isempty(ct2)
                    ct2=ct{i,j}{k1,2};
                else
                    for k=1:length(ct{i,j}{k1,2})
                        if isempty(ct2{k})&&~isempty(ct{i,j}{k1,2}{k})
                            ct2{k}=ct{i,j}{k1,2}{k};
                        end
                        if ~isempty(ct2{k})&&isempty(ct{i,j}{k1,2}{k})
                            ct2{k}=ct2{k}+ct2{k}*0;
                        end  
                        if ~isempty(ct2{k})&&~isempty(ct{i,j}{k1,2}{k})
                            ct2{k}=ct2{k}+ct{i,j}{k1,2}{k};
                            ctt2(k)=ctt2(k)+1;
                        end 
                    end
                end
            end
        end
        
        for k=1:length(ct1)
            ct1{k}=ct1{k}/ctt1(k);
        end
        for k=1:length(ct2)
            ct2{k}=ct2{k}/ctt2(k);
        end
        
        ct_all_dir{1}{i,j}=ct1;
        ct_all_dir{2}{i,j}=ct2;
    end
end

%% ctime
for i=1:size(ctime,1)
    for j=1:size(ctime,2)
        ctime1=[];
        ctime2=[];
        ctimet1=1;
        ctimet2=1;
        for k1=1:size(ctime{i,j},1)
            if ~isempty(ctime{i,j}{k1,1})
                if isempty(ctime1)
                    ctime1=ctime{i,j}{k1,1};
                else
                    ctime1=ctime1+ctime{i,j}{k1,1};
                    ctimet1=ctimet1+1;
                end
            end
            if ~isempty(ctime{i,j}{k1,2})
                if isempty(ctime2)
                    ctime2=ctime{i,j}{k1,2};
                else
                    ctime2=ctime2+ctime{i,j}{k1,2};
                    ctimet2=ctimet2+1;
                end
            end
        end
        
        ctime1=ctime1/ctimet1;
        ctime2=ctime2/ctimet2;
        
        ctime_all_dir{1}{i,j}=ctime1;
        ctime_all_dir{2}{i,j}=ctime2;
    end
end