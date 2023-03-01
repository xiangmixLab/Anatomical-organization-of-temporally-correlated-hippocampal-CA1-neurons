function CM = consensusKmeans_adapted_112722(neuron0,K,N)
%% performing clustering using kmeans+consensus clustering method
% using neuron.trace
rng default  % For reproducibility
CM = zeros(size(neuron0.C,1)); 
options = statset('UseParallel',1);

thresh=3*std(C_to_peakS(neuron0.C),[],2);

for i = 1:N
    datai = neuron0.C;
    datai(datai<thresh)=0;
    
    try
    group = kmeans(datai,K,'Distance','correlation','rep',10,'disp','off','Options',options);
    group = [1:length(group), group];
    [~,idx] = sort(group(:,1)); group = group(idx,2);
    for k = 1:length(unique(group))
        comb = nchoosek(find(group == k),2);
        if length(comb) == 1
            continue;
        end
        linearInd = sub2ind(size(CM), comb(:,1),comb(:,2));
        CM(linearInd) = CM(linearInd) + 1;% if these cells are in same group, similarity +1
    end
    catch
        continue;
    end
end
CM = max(CM,CM')/N;

% cgo = clustergram(CM,'Standardize',3,'Linkage','complete');
% set(cgo,'Colormap',redbluecmap);