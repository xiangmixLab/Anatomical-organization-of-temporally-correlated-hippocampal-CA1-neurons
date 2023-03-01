function [neuron_dir,behav_dir]=laps_avg_neuroDat(neuron_laps,behav_laps)

%% neuron merge
neuron_dir={};
behav_dir={};
for i=1:size(neuron_laps,1)
    for j=1:size(neuron_laps,2)
        nC1=[];
        nC2=[];
        nt1=[];
        nt2=[];
        nC1_ctt=1;
        nC2_ctt=1;

        bp1=[];
        bt1=[];
        bp2=[];
        bt2=[];
        for k1=1:size(neuron_laps{i,j},1)
            % dir1
            if ~isempty(neuron_laps{i,j}{k1,1})
                if isempty(nC1)
                    nC1=neuron_laps{i,j}{k1,1}.C;
                    nt1=neuron_laps{i,j}{k1,1}.time;
                    bp1=behav_laps{i,j}{k1,1}.position;
                    bt1=behav_laps{i,j}{k1,1}.time;
                else
                    nC1=[nC1,neuron_laps{i,j}{k1,1}.C];
                    nt1=[nt1;neuron_laps{i,j}{k1,1}.time];
                    bp1=[bp1;behav_laps{i,j}{k1,1}.position];
                    bt1=[bt1;behav_laps{i,j}{k1,1}.time];
                    nC1_ctt=nC1_ctt+1;
                end
            end
            % dir2
            if ~isempty(neuron_laps{i,j}{k1,2})
                if isempty(nC2)
                    nC2=neuron_laps{i,j}{k1,2}.C;
                    nt2=neuron_laps{i,j}{k1,2}.time;
                    bp2=behav_laps{i,j}{k1,2}.position;
                    bt2=behav_laps{i,j}{k1,2}.time;
                else
                    nC2=[nC2,neuron_laps{i,j}{k1,2}.C];
                    nt2=[nt2;neuron_laps{i,j}{k1,2}.time];
                    bp2=[bp2;behav_laps{i,j}{k1,2}.position];
                    bt2=[bt2;behav_laps{i,j}{k1,2}.time];
                    nC2_ctt=nC2_ctt+1;
                end         
            end
        end
        
        nS1=C_to_peakS(nC1);
        nS2=C_to_peakS(nC1);
        
        neuron_dir{1}{i,j}.C=nC1;
        neuron_dir{1}{i,j}.S=nS1;
        neuron_dir{1}{i,j}.time=nt1;
        neuron_dir{2}{i,j}.C=nC2;
        neuron_dir{2}{i,j}.S=nS2;
        neuron_dir{2}{i,j}.time=nt2;
        
        behav_dir{1}{i,j}.position=bp1;
        behav_dir{1}{i,j}.time=bt1;
        behav_dir{2}{i,j}.position=bp2;
        behav_dir{2}{i,j}.time=bt2;        
    end
end