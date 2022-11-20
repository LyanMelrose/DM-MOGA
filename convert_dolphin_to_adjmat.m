function [ adj_mat,clu_assignment_real ] = convert_dolphin_to_adjmat(network,community)
% Convert the dolphin network data to adjacency matrix

[numedges,temp] = size(network);
clear temp;
numnodes = max(network(:,1))+1;

adj_mat = zeros(numnodes,numnodes);
for i = 1:numedges
    adj_mat(network(i,1)+1,network(i,2)+1) = 1;
end

adj_mat = adj_mat + adj_mat';

clu_assignment_real = community(:,2)';

end