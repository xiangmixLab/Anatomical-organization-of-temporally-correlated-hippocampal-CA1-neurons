% calculate the avg. ratemap correlations across multiple trials with
% pairwise manner

% all_fr: rows--mice, columns--trials, each cell--list of ratemaps of
% multiple neurons

function [place_cells,avg_corr]=pc_by_multi_trial_correlations(all_fr_ori,pc_corr_thresh)
    avg_corr={};
    place_cells={};
    for i=1:size(all_fr_ori)    

        all_fr=all_fr_ori(i,:);
        
        for j=1:length(all_fr)
            for k=1:length(all_fr{j})
                all_fr{j}{k}=filter2DMatrices(all_fr{j}{k},1);
            end
        end

        ratemap_corr_val=[];
        for k=1:length(all_fr{1})
            all_fr_m={};

            for k2=1:size(all_fr,2)
                all_fr_m{k2,1}=reshape(all_fr{k2}{k},1,[]);
%                 all_fr_m{k2,1}(all_fr_m{k2,1}==0)=[];
            end

            all_fr_m=all_fr_m(~cellfun('isempty',all_fr_m));
            
            for k1=2:length(all_fr_m)
                if ~isempty(all_fr_m{k1,1})
                    all_fr_m{k1,1}=resample(all_fr_m{k1,1},length(all_fr_m{1,1}),length(all_fr_m{k1,1}));
                end
            end
            
            all_fr_m1=cell2mat(all_fr_m);

            current_neuron_corr=1-pdist(all_fr_m1,'corr');

            ratemap_corr_val(k,1)=mean(current_neuron_corr);
        end
        avg_corr{i,1}=ratemap_corr_val;
        place_cells{i,1}=find(ratemap_corr_val>pc_corr_thresh);
    end
    