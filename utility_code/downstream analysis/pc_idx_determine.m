function all_pc_inuse=pc_idx_determine(all_pc_dir1,all_neuron_num,tk_range,AD_idx_t,infoscore_type,cond_labels)

pc_idx_dir1={};
for i=tk_range

    pc_idx_dir1{i,1}=zeros(all_neuron_num(i),1);
    pc_idx_dir1{i,2}=zeros(all_neuron_num(i),1);

    t1={};
    t2={};
    
    ctt1=1;
    ctt2=1;
    for k=[1:4]
        t11=all_pc_dir1{i,k}{infoscore_type};

        if ~isempty(findstr(cond_labels{k},'hor'))
            t1{ctt1}=t11; 
            ctt1=ctt1+1;
        end
        if ~isempty(findstr(cond_labels{k},'vet'))
            t2{ctt2}=t11; 
            ctt2=ctt2+1;
        end        
    end
    pc_idx_dir1{i,1}(unique([t1{1};t1{2};t1{3}]))=1;
    pc_idx_dir1{i,2}(unique([t2{1}]))=1;

end

all_pc_inuse={[cell2mat(pc_idx_dir1(AD_idx_t==1,1));cell2mat(pc_idx_dir1(AD_idx_t==1,2))],[cell2mat(pc_idx_dir1(AD_idx_t==2,1));cell2mat(pc_idx_dir1(AD_idx_t==2,2))],[cell2mat(pc_idx_dir1(AD_idx_t==3,1));cell2mat(pc_idx_dir1(AD_idx_t==3,2))],[cell2mat(pc_idx_dir1(AD_idx_t==4,1));cell2mat(pc_idx_dir1(AD_idx_t==4,2))]};

