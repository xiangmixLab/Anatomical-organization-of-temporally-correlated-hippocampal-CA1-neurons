% the direct output of TUnCaT correspond to C_raw of CNMF-E result, we
% deconvolte them to C

function suppl17_CNMFE_style_data_generation(dpath)

    deconv_options = struct('type', 'ar1', ... % model of the calcium traces. {'ar1', 'ar2'}
                            'method', 'thresholded', ... % method for running deconvolution {'foopsi', 'constrained', 'thresholded'} % follow the most original CNMFE BY suoqin, use thresholded instead of foopsi
                            'smin', -5, ...         % minimum spike size. When the value is negative, the actual threshold is abs(smin)*noise level
                            'optimize_pars', true, ...  % optimize AR coefficients
                            'optimize_b', true, ...% optimize the baseline);
                            'max_tau', 100);    % maximum decay time (unit: frame);

    for i=1:length(dpath)
        load([dpath{i},'\','neuronIndividuals_new.mat'])
        load([dpath{i},'\','unmixed_traces','\','alpha= 1.000','\','final_concatenate_stack.mat'])

        % CNMFE sytle decov
        ci=[];
        si=[];
        parfor j=1:size(traces_nmfdemix,2)
            [ci(j,:), si(j,:)] = deconvolveCa(traces_nmfdemix(:,j), deconv_options);
        end

        startIdx=1;
        neuronIndividuals_new_tuncat={};
        for j=1:length(neuronIndividuals_new)
            size_C=size(neuronIndividuals_new{j}.C,2);
            n=neuronIndividuals_new{j}.copy;
            n.C=ci(:,startIdx:startIdx+size_C-1);
            n.S=si(:,startIdx:startIdx+size_C-1);
            n.C_raw=traces_nmfdemix(startIdx:startIdx+size_C-1,:)';
            neuronIndividuals_new_tuncat{j}=n;
            startIdx=startIdx+size_C;
        end

        save([dpath{i},'\','neuronIndividuals_new_tuncat.mat'],'neuronIndividuals_new_tuncat','-v7.3');
    end

