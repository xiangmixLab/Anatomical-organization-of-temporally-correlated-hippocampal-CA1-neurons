behavfile_list={}
for o1=1:length(vname)
    filename1=[vname{o1},'_',behavprefix];
    behavfile_list{o1}=[filename1,'_','Behav.mat'];
end