function [graph_perClust,graph_overall,idx_clust,idx_overall]=clust_graph(neuron,group,corrCut)

ugp=unique(group);

graph_perClust={};
idx_clust={};
for i=1:length(ugp)
    idx=group==ugp(i);
    B1=squareform(1-pdist(neuron.C(idx,:),'correlation'));
    B1(B1<corrCut)=0;
    
    parent={};
    child={};

    idx=[];
    for j1=1:size(B1,1)-1
        for j2=i+1:size(B1,2)
            if B1(j1,j2)>0
                parent=[parent,{num2str(j1)}];
                child=[child,{num2str(j2)}];
                idx=[idx1;j1;j2];
            end
        end
    end

    g1=graph(parent,child);
    graph_perClust{i}=g1;
    idx_clust{i}=unique(idx);
end

% overall graph
B1=squareform(1-pdist(neuron.C,'correlation'));
B1(B1<corrCut)=0;

parent={};
child={};

idx=[];
for j1=1:size(B1,1)-1
    for j2=i+1:size(B1,2)
        if B1(j1,j2)>0
            parent=[parent,{num2str(j1)}];
            child=[child,{num2str(j2)}];
            idx=[idx1;j1;j2];
        end
    end
end

g1=graph(parent1,child1);
graph_overall=g1;
idx_overall=unique(idx);
