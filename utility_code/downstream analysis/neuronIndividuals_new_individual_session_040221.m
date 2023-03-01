%% HDAC AD processing script

function neuronIndividuals_new_individual_session_040221(foldernamet,timestampname,mscamid,behavcamid,mrange)

    %% start
    t1=tic;
    for ikk=mrange  

        cd(foldernamet{ikk})
        
        % read data
        load('all_neurons.mat');

        % split neuron data
        mscamid1=mscamid{ikk};
        behavcamid1=behavcamid{ikk};
        neuronIndividuals_new={};
        for i = 1:length(neurons)
            neuronIndividuals_new{i}=neurons{i}.copy;
            if ~isempty(timestampname{1,i})
                fid=fopen(timestampname{1,i},'r');
                mscamindex=mscamid1(1,i);
                behavindex=behavcamid1(1,i);
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
        save('neuronIndividuals_new_trial.mat','neuronIndividuals_new','-v7.3');
        disp(['fin ',num2str(ikk),' use ',num2str(toc(t1)),' sec']);
    end