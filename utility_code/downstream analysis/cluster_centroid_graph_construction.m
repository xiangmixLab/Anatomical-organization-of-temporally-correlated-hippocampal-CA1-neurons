function G1=centroid_graph_construction(cen)


distt=squareform(pdist(cen));
distt(distt==0)=nan;
distt_neighbor_thresh=35; % only count the neurons within 35um

ctt=1;
node_list=[];
for k1=1:size(cen,1)
    
    dir_cart=cen-cen(k1,:);
    [theta,rho] = cart2pol(dir_cart(:,1),dir_cart(:,2)); % convert to polar coordination sys
    theta=rad2deg(theta);

    for k2=1:length(rho) % for each neighborhood point (include itself)
        non_neighbor_sign=0;
        for k21=[min(1:k2-1,1),min(k2+1:length(rho),length(rho))] % check all other points
            rho_curr=rho(k2);
            rho_other=rho(k21);

            % determine if cen(k1,:) and cen(k2,:) are neighbors:
            % 1: judge if there's a point between current point pair 
            % that is, cen(k1,:) and cen(k2,:), check if cen(k21,:) is nearer than cen(k2,:) 
            % and whether cen(k21,:) is within +/-10 degree of the path between cen(k1,:) and cen(k2,:)
            % 2: distance between cen(k1,:) and cen(k2,:) smaller than 5%
            % of lower distances between cen(k1,:) and others
            % 3: cen(k2,:) is not cen(k1,:)
            if (rho_other<rho_curr && abs(theta(k2)-theta(k21))<=10) || rho_curr>distt_neighbor_thresh || rho_curr==0
                non_neighbor_sign=1; % if so, cen(k1,:) and cen(k2,:) are not neighbors
                break;
            end
        end
        if non_neighbor_sign==0 
            node_list(ctt,1)=k1;
            node_list(ctt,2)=k2; 
            ctt=ctt+1;
        end
        % if no points in between cen(k1,:) and cen(k2,:), then this two points are neighboors. 
    end
end

dup_idx=zeros(size(node_list,1),1);
for i=1:size(node_list,1)-1
    for j=i+1:size(node_list,1)
        if (node_list(i,1)==node_list(j,1)&&node_list(i,2)==node_list(j,2))||(node_list(i,1)==node_list(j,2)&&node_list(i,2)==node_list(j,1))
            dup_idx(j,1)=1;
        end
    end
end

node_list(dup_idx==1,:)=[];

G1=graph(node_list(:,1),node_list(:,2)); % in this way we create a graph in which each node only connect to its neighboors

% p=plot(G1);
% p.MarkerSize=5;
% p.XData=cen(:,1);
% p.YData=cen(:,2);
