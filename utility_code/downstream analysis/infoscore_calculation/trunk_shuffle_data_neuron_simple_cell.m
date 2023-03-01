function neuron0=trunk_shuffle_data_neuron_simple_cell(C,S,time,cell_idx)

C1=trunk_shuffle_data_pre_split(C);
S1=trunk_shuffle_data_pre_split(S);
neuron0.C=C1(cell_idx,:);
neuron0.S=S1(cell_idx,:);
neuron0.time=time;
