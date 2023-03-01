function all_neuron1=trunk_shuffle_data_laps(all_neuron,trunk_num)

all_neuron1={};
for k1=1:size(all_neuron,1)
    for k2=1:size(all_neuron,2)
        
        if ~isempty(all_neuron{k1,k2})
            dat1=all_neuron{k1,k2}.C;
            dat2=all_neuron{k1,k2}.C_raw;
            dat3=all_neuron{k1,k2}.S;

            interval=round(size(dat1,2)/(trunk_num));
            interval_list=[ones(1,floor(size(dat1,2)/interval))*interval,size(dat1,2)-floor(size(dat1,2)/interval)*interval];
            interval_list(interval_list==0)=[]; % the last one is 0 inperfact case

            dat_trunk1=mat2cell(dat1',interval_list); %mat2cell(dat,[row_sizes_of_the_mat_in_current_splitted_cell],[col_sizes_of_the_mat_in_current_splitted_cell])
            dat_trunk2=mat2cell(dat2',interval_list); %mat2cell(dat,[row_sizes_of_the_mat_in_current_splitted_cell],[col_sizes_of_the_mat_in_current_splitted_cell])
            dat_trunk3=mat2cell(dat3',interval_list); %mat2cell(dat,[row_sizes_of_the_mat_in_current_splitted_cell],[col_sizes_of_the_mat_in_current_splitted_cell])

            rand_order=randperm(size(dat_trunk1,1));
            dat_shuf1=cell2mat(dat_trunk1(rand_order));
            dat_shuf2=cell2mat(dat_trunk2(rand_order));
            dat_shuf3=cell2mat(dat_trunk3(rand_order));

            all_neuron1{k1,k2}.C=dat_shuf1';
            all_neuron1{k1,k2}.C_raw=dat_shuf2';
            all_neuron1{k1,k2}.S=dat_shuf3';
            all_neuron1{k1,k2}.time=all_neuron{k1,k2}.time;
        end
    end
end