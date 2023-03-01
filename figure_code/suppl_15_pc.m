% Fig2 place cells

% place cell calculation
foldername={
    'D:\Remapping_linear_track_053019_053119\result_merge\M3411'
    'D:\Remapping_linear_track_053019_053119\result_merge\M3412'
    'D:\Remapping_linear_track_053019_053119\result_merge\M3421F'
    'D:\Remapping_linear_track_053019_053119\result_merge\M3422F'
    'D:\Remapping_linear_track_053019_053119\result_merge\M3424F'
    'D:\Remapping_linear_track_053019_053119\result_merge\M3425F'
    }


behavname={
    'Horizontal1_Behav.mat';
    'Horizontal2_Behav.mat'; 
    'Vertical_Behav.mat'; 
    }

all_behav={};
for i=1:length(foldername)
    for j=1:length(behavname)
        load([foldername{i},'\',behavname{j}]);
        all_behav{i,j}=behav;
        all_behav{i,j}.VidObj=[];
    end
end

% PC
all_pc=cell(6,3);
all_infoscore=cell(6,3);
all_infoscore_norm=cell(6,3);
all_coherence=cell(6,3);
tic;
for i=1:6
    load([foldername{i},'\','neuronIndividuals_new.mat']);
    for j=1:3
    
        [place_cells,infoScore,infoScore_norm,coherencee] = permutingSpike_adapt_040821(neuronIndividuals_new{j},all_behav{i,j}.position,all_behav{i,j}.time,'S',0,10,5,'all',0.4);  

        all_pc{i,j}=place_cells;
        all_infoscore{i,j}=infoScore;
        all_infoscore_norm{i,j}=infoScore_norm;
        all_coherence{i,j}=coherencee;

    end
end
toc;

% pc by laps
all_pc_laps=cell(6,3);
all_infoscore_laps=cell(6,3);
all_infoscore_norm_laps=cell(6,3);
all_coherence_laps=cell(6,3);
tic;
for i=1:6
    for j=1:3
        for k1=1:size(nnew{i,j},1)
            for k2=1:size(nnew{i,j},2)
                if ~isempty(nnew{i,j}{k1,k2})
                    [place_cells,infoScore,infoScore_norm,coherencee] = permutingSpike_adapt_040821(nnew{i,j}{k1,k2},behavnew{i,j}{k1,k2}.position,behavnew{i,j}{k1,k2}.time,'S',0,10,1,'all',0.4);  

                    all_pc_laps{i,j}{k1,k2}=place_cells;
                    all_infoscore_laps{i,j}{k1,k2}=infoScore;
                    all_infoscore_norm_laps{i,j}{k1,k2}=infoScore_norm;
                    all_coherence_laps{i,j}{k1,k2}=coherencee;
                end
            end
        end

    end
end
toc;

% pc by laps, ratemap avg
all_pc_laps_avg=cell(6,3);
all_infoscore_laps_avg=cell(6,3);
all_infoscore_norm_laps_avg=cell(6,3);
all_coherence_laps_avg=cell(6,3);
tic;

for i=1:6
    for j=1:3

        [place_cells,infoScore,infoScore_norm,coherencee] = permutingSpike_adapt_laps_040221(nnew{i,j},behavnew{i,j},'S',0,15,nnum(i,j),1,'all',0.4);  

        all_pc_laps_avg{i,j}=place_cells;
        all_infoscore_laps_avg{i,j}=infoScore;
        all_infoscore_norm_laps_avg{i,j}=infoScore_norm;
        all_coherence_laps_avg{i,j}=coherencee;


    end
end

% pc, original, zscore ratemap
all_pc_z=cell(6,3);
all_infoscore_z=cell(6,3);
all_infoscore_norm_z=cell(6,3);
all_coherence_z=cell(6,3);
tic;
for i=1:6
    load([foldername{i},'\','neuronIndividuals_new.mat']);
    for j=1:3
    
        [place_cells,infoScore,infoScore_norm,coherencee] = permutingSpike_adapt_041222_LT(neuronIndividuals_new{j},all_behav{i,j}.position,all_behav{i,j}.time,'S',0,10,5,'all',0.4);  

        all_pc_z{i,j}=place_cells;
        all_infoscore_z{i,j}=infoScore;
        all_infoscore_norm_z{i,j}=infoScore_norm;
        all_coherence_z{i,j}=coherencee;

    end
end

