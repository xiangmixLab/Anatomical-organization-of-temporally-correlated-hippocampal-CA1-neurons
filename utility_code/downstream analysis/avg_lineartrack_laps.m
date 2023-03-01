function [fr_all_laps_avg,fr_all_laps_avg_line,ctime_all_laps_avg]=avg_lineartrack_laps(fr_all_laps,ctime_all_laps,all_neuron_num)

fr_all_laps_avg={};
fr_all_laps_avg_line={};
ctime_all_laps_avg={};

for tk=1:size(fr_all_laps,1)
    for i=1:size(fr_all_laps,2)
        fr_all_laps_temp1=fr_all_laps{tk,i}(:,1);
        fr_all_laps_temp2=fr_all_laps{tk,i}(:,2);
        
        % max size determine
        size_mat=[];
        for j=1:size(fr_all_laps_temp1,1)
            if ~isempty(fr_all_laps_temp1{j,1})
                for j1=1:size(fr_all_laps_temp1{j,1},2)
                    size_mat=[size_mat;size(fr_all_laps_temp1{j,1}{j1})];
                end
            end
            if ~isempty(fr_all_laps_temp2{j,1})
                for j1=1:size(fr_all_laps_temp2{j,1},2)
                    size_mat=[size_mat;size(fr_all_laps_temp2{j,1}{j1})];
                end
            end
        end
        maxsize=max(size_mat,[],1);
        
        % reshape and fill in 0
        for j=1:size(fr_all_laps_temp1,1)
            if ~isempty(fr_all_laps_temp1{j})
                for j1=1:size(fr_all_laps_temp1{j},2)
                    if ~isempty(fr_all_laps_temp1{j}{j1})
                        fr_all_laps_temp1{j}{j1}=imresize(fr_all_laps_temp1{j}{j1},maxsize);
                    else
                        fr_all_laps_temp1{j}{j1}=zeros(maxsize(1),maxsize(2));
                    end
                end
            end
            if ~isempty(fr_all_laps_temp2{j})
                for j1=1:size(fr_all_laps_temp2{j},2)
                    if ~isempty(fr_all_laps_temp2{j}{j1})
                        fr_all_laps_temp2{j}{j1}=imresize(fr_all_laps_temp2{j}{j1},maxsize);
                    else
                        fr_all_laps_temp2{j}{j1}=zeros(maxsize(1),maxsize(2));
                    end
                end
            end
        end
        
        % average
        fr_avg_t1={};
        fr_avg_t2={};    
        fr_avg_l_t1=[];
        fr_avg_l_t2=[];
        for j1=1:all_neuron_num(tk)
            temp_avg_rm1=zeros(maxsize(1),maxsize(2));
            temp_avg_rm2=zeros(maxsize(1),maxsize(2));
            for j=1:size(fr_all_laps_temp1,1)
                if ~isempty(fr_all_laps_temp1{j})
                     temp_avg_rm1= temp_avg_rm1+fr_all_laps_temp1{j}{j1};
                end
                if ~isempty(fr_all_laps_temp2{j})
                     temp_avg_rm2= temp_avg_rm2+fr_all_laps_temp2{j}{j1};     
                end
            end
            fr_avg_t1{j1,1}=temp_avg_rm1/(size(fr_all_laps_temp1,1)/2);
            fr_avg_t2{j1,1}=temp_avg_rm2/(size(fr_all_laps_temp2,1)/2);
            
            if size(fr_avg_t1{j1,1},1)<size(fr_avg_t1{j1,1},2) % horizontal track
                fr_avg_l_t1(j1,:)=nansum(fr_avg_t1{j1,1},1);
                fr_avg_l_t2(j1,:)=nansum(fr_avg_t2{j1,1},1);
            else % vertical track
                fr_avg_l_t1(j1,:)=nansum(fr_avg_t1{j1,1},2)'; % 90 degree dir rotation
                fr_avg_l_t2(j1,:)=nansum(fr_avg_t2{j1,1},2)';  
            end
        end
        
        temp_avg_ctime1=zeros(maxsize(1),maxsize(2));
        temp_avg_ctime2=zeros(maxsize(1),maxsize(2));
        for j=1:size(ctime_all_laps{tk,i},1)
            if ~isempty(ctime_all_laps{tk,i}{j,1})
                 temp_avg_rm1= temp_avg_rm1+ctime_all_laps{tk,i}{j,1};
            end
            if ~isempty(ctime_all_laps{tk,i}{j,2})
                 temp_avg_rm2= temp_avg_rm2+ctime_all_laps{tk,i}{j,2};   
            end
        end
        ctime_avg_t1=temp_avg_rm1/(size(fr_all_laps_temp1,1)/2);
        ctime_avg_t2=temp_avg_rm2/(size(fr_all_laps_temp2,1)/2);
        
        fr_all_laps_avg{tk,i}{1}=fr_avg_t1;
        fr_all_laps_avg{tk,i}{2}=fr_avg_t2;
        fr_all_laps_avg_line{tk,i}{1}=fr_avg_l_t1;
        fr_all_laps_avg_line{tk,i}{2}=fr_avg_l_t2;
        ctime_all_laps_avg{tk,i}{1}=ctime_avg_t1;
        ctime_all_laps_avg{tk,i}{2}=ctime_avg_t2;
    end
end