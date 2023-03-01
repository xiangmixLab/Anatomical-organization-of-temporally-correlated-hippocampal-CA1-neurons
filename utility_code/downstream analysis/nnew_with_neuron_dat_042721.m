%% HDAC AD processing script

function neuronIndividuals_new=nnew_with_neuron_dat_042721(neurons,timestampname,mscamid,behavcamid)

    %% start
    t1=tic;

    % split neuron data
    mscamid1=mscamid;
    behavcamid1=behavcamid;
    neuronIndividuals_new={};
    for i = 1:length(neurons)
        neuronIndividuals_new{i}=neurons{i}.copy;
        if ~isempty(timestampname{1,i})
            fid=fopen(timestampname{1,i},'r');
            mscamindex=mscamid1;
            behavindex=behavcamid1;
            timedata = textscan(fid, '%d%d%d%d', 'HeaderLines', 1);
            t = timedata{3}(timedata{1} == mscamindex);t(1) = 0;   %% make sure the miniscope is for USB port 1 %%behav cam port 0
            %%othewise %%t = timedata{3}(timedata{1} == 0);t(1) = 0;
            idxt = find(diff(t)<=0);
            t(idxt+1) = t(idxt)+1;
            % the first is trainning and the second is testing
            neuronIndividuals_new{i}.time = double(t);
        end
    end

    % save
    disp(['fin ',' use ',num2str(toc(t1)),' sec']);
