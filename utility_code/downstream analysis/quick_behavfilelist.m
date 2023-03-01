behavfile_list={};
behavVidprefix='behavCam1_';
for i=1:length(orilocation)
    behavfile_list{i}=[vname{i},'_',behavVidprefix];
    behavfile_list{i}=[behavfile_list{i},'_','Behav.mat'];
end