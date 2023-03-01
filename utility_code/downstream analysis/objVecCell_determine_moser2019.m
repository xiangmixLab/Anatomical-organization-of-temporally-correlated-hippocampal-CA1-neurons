function [objVecCell_idx,robj1,robj1_shuf]=objVecCell_determine_moser2019(fr1,fr2,obj_loc1,obj_loc2)

% fr1: ratemap at first obj loc
% fr2: ratemap at second obj loc
dis_obj_graph1={};
dis_obj_graph2={};
for k=1:length(fr1)
    dis_obj_graph1{k,1}=obj_dis_ang_plot_convert(fr1{k,1},obj_loc1);
    dis_obj_graph2{k,1}=obj_dis_ang_plot_convert(fr2{k,1},obj_loc2);
end

robj1=rateMap_correlation(dis_obj_graph1,dis_obj_graph2,ones(size(dis_obj_graph1{1,1})),ones(size(dis_obj_graph2{1,1})),1);

robj1_shuf=[];

parfor i=1:100
    dshuf1=dis_obj_graph1;
    dshuf2=dis_obj_graph2;
    for k=1:length(fr1)
        dshuf1{k,1}=reshape(dshuf1{k,1}(randperm(size(dshuf1{k,1},1)*size(dshuf1{k,1},2))),size(dshuf1{k,1},1),size(dshuf1{k,1},2));
        dshuf2{k,1}=reshape(dshuf2{k,1}(randperm(size(dshuf2{k,1},1)*size(dshuf2{k,1},2))),size(dshuf2{k,1},1),size(dshuf2{k,1},2));
    end

    robj1_shuf(:,i)=rateMap_correlation(dshuf1,dshuf2,ones(size(dshuf1{1,1})),ones(size(dshuf2{1,1})),1);

end

robj1_thresh=quantile(robj1_shuf,0.95,2);
objVecCell_idx=(robj1'>robj1_thresh);
objVecCell_idx=(robj1'>robj1_thresh).*(robj1'>=0.5);