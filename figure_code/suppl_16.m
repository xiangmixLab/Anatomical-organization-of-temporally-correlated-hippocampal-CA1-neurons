% number of fields of cells
foldername_all={
    {
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_1_3_tri_cir_sqr\M3411'	
    'D:\\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_1_3_tri_cir_sqr\M3412'	
    'D:\\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_1_3_tri_cir_sqr\M3421F'	
    'D:\\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_1_3_tri_cir_sqr\M3422F'	
    'D:\\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_1_3_tri_cir_sqr\M3424F'	
    'D:\\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_1_3_tri_cir_sqr\M3425F'	
    };
    {
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_2_cir_rec\M3411'		
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_2_cir_rec\M3412'		
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_2_cir_rec\M3413'		
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_2_cir_rec\M3414'		
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_2_cir_rec\M3415'		
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_2_cir_rec\M3421F'		
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_2_cir_rec\M3422F'		
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_2_cir_rec\M3423F'		
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_2_cir_rec\M3424F'		
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_2_cir_rec\M3425F'		
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_2_cir_rec\M3426F'		
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_2_cir_rec\M3427F'		
    };

    {
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_3_barrier\101F'
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_3_barrier\102F'
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_3_barrier\103F'
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_3_barrier\840F'
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_3_barrier\2833M'
    };

    {
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_3_linear_track\M3411'
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_3_linear_track\M3412'
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_3_linear_track\M3421F'
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_3_linear_track\M3422F'
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_3_linear_track\M3424F'
    'D:\Xu_clusterting_paper_prep11_2020\final data\final_neuron_behav_data\Fig_3_linear_track\M3425F'
    }
    }

behavName_all={
    {
    'triangle1_Behav.mat';
    'circle1_Behav.mat';
    'square1_Behav.mat';   
    'circle2_Behav.mat';
    'square2_Behav.mat';   
    'triangle2_Behav.mat';
    };

    {
    'square_Behav.mat';
    'circle_Behav.mat';
    };

    {
    'training_precue_Behav.mat'		
    'training_cue_Behav.mat'		
    'training_postcue_Behav.mat'	
    'barrier_precue_Behav.mat'
    'barrier_cue_Behav.mat'		
    'barrier_postcue_Behav.mat'		
    };

    {
    'Horizontal1_Behav.mat';
    'Horizontal2_Behav.mat'; 
    'Vertical_Behav.mat'; 
    }
    }

cond={
    [1:6];
    [1:2];
    [1:6];
    [1:3];
    }
% all field nums
all_field_nums={};
firingrateAll={};
countTime={};
for i=1:length(foldername_all)
    for j=1:length(foldername_all{i})
        load([foldername_all{i}{j},'\','neuronIndividuals_new.mat'])
        for k=cond{i}
            load([foldername_all{i}{j},'\',behavName_all{i}{k}]);
            
            neuron=neuronIndividuals_new{k}.copy;
            behavpos=behav.position;
            behavtime=behav.time;
            maxbehavROI=behav.ROI;
            binsize=10;
            thresh=3*std(C_to_peakS(neuron.C),[],2);
            countTimeThresh=[0 inf];
            small_velo=10;

            [firingrateAll{i}{j,k},countAll,~,countTime{i}{j,k},~,~,bininfo,cpt] = calculatingCellSpatialForSingleData_040321(neuron,behavpos,behavtime,maxbehavROI,binsize,1:size(neuron.C,1),thresh,'S',[],[],countTimeThresh,small_velo);
            
            numFields=field_number_calculation(firingrateAll{i}{j,k});

            all_field_nums{i}{j,k}=numFields';
        end
    end
end

%% illustration
dat1=[cell2mat(all_field_nums{1}(:,1));cell2mat(all_field_nums{1}(:,2));cell2mat(all_field_nums{1}(:,3));cell2mat(all_field_nums{1}(:,4));cell2mat(all_field_nums{1}(:,5));cell2mat(all_field_nums{1}(:,6))];
dat2=[cell2mat(all_field_nums{2}(:,1));cell2mat(all_field_nums{2}(:,2))];
dat3=[cell2mat(all_field_nums{3}(:,1));cell2mat(all_field_nums{3}(:,2));cell2mat(all_field_nums{3}(:,3));cell2mat(all_field_nums{3}(:,4));cell2mat(all_field_nums{3}(:,5));cell2mat(all_field_nums{3}(:,6))];
dat4=[cell2mat(all_field_nums{4}(:,1));cell2mat(all_field_nums{4}(:,2));cell2mat(all_field_nums{4}(:,3))];

dat1(dat1==0)=[];
dat2(dat2==0)=[];
dat3(dat3==0)=[];
dat4(dat4==0)=[];

[dat11,e]=histcounts(dat1,'binWidth',1);
[dat21,e]=histcounts(dat2,'binWidth',1);
[dat31,e]=histcounts(dat3,'binWidth',1);
[dat41,e]=histcounts(dat4,'binWidth',1);

subplot(141)
bar(dat11);hold on;
disp(num2str(dat11./sum(dat11)))

subplot(142)
bar(dat21);hold on;
disp(num2str(dat21./sum(dat21)))

subplot(143)
bar(dat31);hold on;
disp(num2str(dat31./sum(dat31)))

subplot(144)
bar(dat41);hold on;
disp(num2str(dat41./sum(dat41)))

%% illustration: cells with different num of fields
all_dat={dat1,dat2,dat3,dat4}
for i=1:4
    totalNum=max(all_dat{i});
    figure;
    for p=1:totalNum
        
        idx=[];
        for k1=1:size(all_field_nums{i},1)
            for k2=1:size(all_field_nums{i},2)
                if ~isempty(find(all_field_nums{i}{k1,k2}==p))
                    idx1=find(all_field_nums{i}{k1,k2}==p);
                    idx=[k1,k2,idx1(1)];
                end
            end
        end
        
        subplot(ceil(totalNum/2),2,p);
        [bx2,by2]=boundary_from_binImg(firingrateAll{i}{idx(1),idx(2)}{idx(3)},0.5,1);
        ratemap_plot(firingrateAll{i}{idx(1),idx(2)}{idx(3)},countTime{i}{idx(1),idx(2)},1,0,[]);hold on; 
        for l=1:length(bx2)
            plot(bx2{l},by2{l},'--','color','k'); 
        end
    end
    set(gcf,'renderer','painters');
end

        
