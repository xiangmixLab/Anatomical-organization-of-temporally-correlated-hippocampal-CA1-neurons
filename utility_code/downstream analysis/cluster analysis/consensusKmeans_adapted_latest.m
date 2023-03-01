%% performing clustering using kmeans+consensus clustering method

%% input: 
% neuron0: neuron data
% K: target cluster number
% N: number of bootstrap repeat
%% output CM: consensus matrix representing the times two neurons are assigned to the same clust

function CM = consensusKmeans_adapted_latest(neuron0,K,N)

rng default  % For reproducibility
CM = zeros(size(neuron0.C,1)); 
options = statset('UseParallel',0); % parallel will be done with outer loops
groups={};

for i = 1:N

    oriLeng=size(neuron0.C,2);

    idxx=[1:oriLeng];
    idxx1=randperm(oriLeng);
    idxx(idxx1(1:round(oriLeng*0.5)))=0; % randomly remove 50% of data
    datai = neuron0.C(:,logical(idxx));
    
    datai = interp1(find(idxx~=0),datai',[1:length(idxx)])';
    datai(isnan(datai))=neuron0.C(isnan(datai));

    idxD = find(std(datai,[],2) < 0.01); % identify resampled neurons that have very weak response (very rare by now), make them the same cluster, they will not participate kcc as the function will turn errors
    idxC = setdiff(1:size(datai,1),idxD); 
    datai(idxD,:) = [];
    groupD = [idxD,(K+1)*ones(length(idxD),1)];

    try
        group = kmeans(datai,K,'Distance','correlation','rep',10,'disp','off','Options',options);
        groupC = [idxC(:),group];
        group = [groupC; groupD];
        [~,idx] = sort(group(:,1)); 
        group = group(idx,2);
        groups{i}=group;
    catch
        disp([num2str(K),', round',num2str(i),' ','error']);
        continue;
    end
end

for i=1:N
    for k = 1:length(unique(groups{i}))
            comb = nchoosek(find(groups{i} == k),2);
            if length(comb) == 1
                continue;
            end
            linearInd = sub2ind(size(CM), comb(:,1),comb(:,2));
            CM(linearInd) = CM(linearInd) + 1;% if these cells are in same group, similarity +1
    end
end
CM = max(CM,CM')/N;
