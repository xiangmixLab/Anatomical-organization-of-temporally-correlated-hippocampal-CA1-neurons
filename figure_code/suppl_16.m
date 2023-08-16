% number of fields of cells
foldername_all={
    foldername_multiGeo;
    foldername_fig2;
    foldername_AI163;
    foldername_LT;
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
            load([foldername_all{i}{j},'\','behav.mat']);
            
            neuron=neuronIndividuals_new{k};
            behav=behavIndividuals{k};
            behavpos=behav.position;
            behavtime=behav.time;
            binsize=10;
            thresh=3*std(C_to_peakS(neuron.C),[],2);
            countTimeThresh=[0 inf];
            small_velo=10;

            [firingrateAll{i}{j,k},countAll,~,countTime{i}{j,k},~,~,bininfo,cpt] = calculatingCellSpatialForSingleData_040321(neuron,behavpos,behavtime,[0,0,max(behavpos,[],1)],binsize,1:size(neuron.C,1),thresh,'S',[],[],countTimeThresh,small_velo);
            
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

