function suppl17_mask_gen(dpath)

for i=1:length(dpath)
    load([dpath{i},'\','neuronIndividuals_new.mat']);
    FinalMasks=stacked_footprint_gen(neuronIndividuals_new{1});
    save([dpath{i},'\','FinalMasks_',fn,'.mat'],'FinalMasks')
end