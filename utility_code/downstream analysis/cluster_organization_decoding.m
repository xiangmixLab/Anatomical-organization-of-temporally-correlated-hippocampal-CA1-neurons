function cluster_organization_decoding(behav,neuron,gp_select,neuroDat)

if ~isempty(neuroDat)
    % 50% trial split
    [all_neuron_decode_tr,all_neuron_decode_ts,all_behav_decode_tr,all_behav_decode_ts]=decoding_split_trial(neuron,behav,Fs,[]);

    %% 1: for training: get bin time map and corresponding cluster organization

    twin=750; %50sec

    group_acrossTime={};
    all_nperiod_acrossTime={};
    all_behav_acrossTime={};

    tic;
    ctt=1;
    for k=1:Fs:size(all_neuron_decode_tr.C,2)-twin % every 1 sec
        n_period=all_neuron_decode_tr;
        n_period.C=n_period.C(:,k:k+twin-1);
        n_period.S=C_to_peakS(n_period.C);

        group_acrossTime{1,ctt}=cluster_determine_by_suoqin_NMF_firstPeakCoph(n_period,100,10,gp_select);
        all_nperiod_acrossTime{1,ctt}=n_period;   
        ctt=ctt+1;
    end
    toc;

    ctt=1;
    ntime=all_neuron_decode_tr.time;
    ntime=interp1(ntime,ntime,all_behav_decode_tr.time);
    for k1=1:75*2:size(all_behav_decode_tr.position,1)-twin*2 % every sec
        all_nperiod_acrossTime{1,ctt}.time=ntime(k1:k1+twin*2-1);
        all_behav_acrossTime{1,ctt}.position=all_behav_decode_tr.position(k1:k1+twin*2-1,:);
        all_behav_acrossTime{1,ctt}.time=all_behav{1,j1}.time(k1:k1+twin*2-1);
        all_behav_acrossTime{1,ctt}.ROI=all_behav{1,j1}.ROI;
        ctt=ctt+1;
    end
end
    


