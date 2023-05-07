% suppl 14: dendrite extraction

%% part 1: generate data
% no dendrite: with_dendrite = false

foldernamestruct=destination_no_dendrite;

for i=[1]
    fd=fileparts(foldernamestruct{i});
    cd(fd);
    tic;
        SCOUT_pipeline_single_031422(foldernamestruct{i},[0],[240,376],false,9,20); % changed to BatchEndoscopeauto_adapted_new_kevin 040420 
%         SCOUT_pipeline_single(foldernamestruct{i},[0],[240,376],false); % changed to BatchEndoscopeauto_adapted_new_kevin 040420 

    toc;
    close all;
end

% dendrite: with_dendrite = true

foldernamestruct=destination_dendrite;

for i=[1:length(foldernamestruct)]
    fd=fileparts(foldernamestruct{i});
    cd(fd);
    tic;
        SCOUT_pipeline_single_031422(foldernamestruct{i},[0],[240,376],true,9,20); % changed to BatchEndoscopeauto_adapted_new_kevin 040420 
    toc;
end

%% handle overlap
foldernamestruct=destination_no_dendrite;

for i=[1:length(foldernamestruct)]
    fd=fileparts(foldernamestruct{i});
    cd(fd);
    tic;
        
    pth=fileparts(foldernamestruct{i});
    load([pth,'\','further_processed_neuron_extraction_final_result.mat']);
%     [neuron0,indices_rem]=overlap_handling(neuron);
    
    [Asize1,~,neuron1]=neuronAreaAndTrim(neuron,400);
    
    neuron=neuron1.copy;
    save([pth,'\','further_processed_neuron_extraction_final_result_sz.mat'],'neuron','-v7.3');

    toc;
    close all;
end

% dendrite: with_dendrite = true

foldernamestruct=destination_dendrite;

for i=[1:length(foldernamestruct)]
    fd=fileparts(foldernamestruct{i});
    cd(fd);
    tic;
    
    pth=fileparts(foldernamestruct{i});
    load([pth,'\','further_processed_neuron_extraction_final_result.mat']);
%     [neuron0,indices_rem]=overlap_handling(neuron);
%     neuron=neuron0.copy;
    [Asize1,~,neuron1]=neuronAreaAndTrim(neuron,400);

    save([pth,'\','further_processed_neuron_extraction_final_result_sz.mat'],'neuron','-v7.3');

    toc;
end

%% part 2: illustration

midx={[6,4,5],[10,13,12],[20,21,23]};

% dendrite

for tk=1:length(midx)
    midxTmp=midx{tk};
    figure;
    for i=1:length(midxTmp)
        subplot(2,length(midxTmp),i);
        pth=fileparts(destination_dendrite{midxTmp(i),1});
        vidName=[pth,'\','vid.mat'];
        neuronName=[pth,'\','further_processed_neuron_extraction_final_result.mat'];

        plot_contour_Proj(vidName,neuronName)

        subplot(2,length(midxTmp),i+length(midxTmp));
        load(neuronName);
        footprint=spatial_footprint_calculation(neuron,0.5);
        imagesc(footprint);
    end
    set(gcf,'renderer','painters')
end


% no dendrite
for tk=1:length(midx)
    midxTmp=midx{tk};
    figure;
    for i=1:length(midxTmp)
        subplot(2,length(midxTmp),i);
        pth=fileparts(destination_no_dendrite{midxTmp(i),1});
        vidName=[pth,'\','vid.mat'];
        neuronName=[pth,'\','further_processed_neuron_extraction_final_result.mat'];

        plot_contour_Proj(vidName,neuronName)

        subplot(2,length(midxTmp),i+length(midxTmp));
        load(neuronName);
        footprint=spatial_footprint_calculation(neuron,0.5);
        imagesc(footprint);
    end
    set(gcf,'renderer','painters')
end

%% imshowpair
figure;
for tk=1:length(midx)
    midxTmp=midx{tk};
    figure;
    for i=1:length(midxTmp)
   
        pth=fileparts(destination_dendrite{midxTmp(i),1});
        neuronName=[pth,'\','further_processed_neuron_extraction_final_result.mat'];
        load(neuronName);
        footprint1=spatial_footprint_calculation(neuron,0.5);
        
        pth=fileparts(destination_no_dendrite{midxTmp(i),1});
        neuronName=[pth,'\','further_processed_neuron_extraction_final_result.mat'];
        load(neuronName);
        footprint2=spatial_footprint_calculation(neuron,0.5);
        
        subplot(1,length(midxTmp),i);
        imshowpair(footprint1,footprint2);
    end
end

%% DIFFERENCE : directly substract footprint, see how much left

for tk=1:length(midx)
    midxTmp=midx{tk};
    figure;
    for i=1:length(midxTmp)
   
%         [diff_footprint,diff_coor,diff_neuron]=footprintDiff(destination_dendrite{midxTmp(i),1},destination_no_dendrite{midxTmp(i),1});
        
%         subplot(2,length(midxTmp),i);
%         imagesc(diff_neuron.Cn); hold on
%         plot_contour_Proj2(diff_neuron);
%         
%         subplot(2,length(midxTmp),i+length(midxTmp));
%         imagesc(diff_footprint)

        pth=fileparts(destination_dendrite{midxTmp(i),1});
        neuronName=[pth,'\','further_processed_neuron_extraction_final_result.mat'];
        load(neuronName);
        footprint1=spatial_footprint_calculation(neuron,0.5)>0;
        
        pth=fileparts(destination_no_dendrite{midxTmp(i),1});
        neuronName=[pth,'\','further_processed_neuron_extraction_final_result.mat'];
        load(neuronName);
        footprint2=spatial_footprint_calculation(neuron,0.5)>0;
        
        diffFP=footprint1-footprint2;
        diffFP(diffFP>0)=2;
        diffFP(diffFP<0)=0;
        
        diff_mask=diffFP~=0;
        diff_mask=bwareaopen(diff_mask,9);
        diffFP=diffFP.*double(diff_mask);
        subplot(1,length(midxTmp),i);
        imagesc(diffFP)
    end
end
