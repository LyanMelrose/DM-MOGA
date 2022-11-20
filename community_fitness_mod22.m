function [ cf ] = community_fitness_mod2( adj_mat,clu_assignment,ll )
% Compute the community fittness for each partition.
% adj_mat: the adjacency matrix of the network.
% clu_assignment: the cluster label vector.

a = 1;  % a parameter in the formula of community fittness.
cf = 0;
ec=0;
clu_num = max(clu_assignment);
for i = 1:clu_num
    s_index = find(clu_assignment == i);
    s = adj_mat(s_index,s_index);
    s_cardinality = length(s_index);
    kins_sum = 0;
    kouts_sum = 0;
    for j = 1:s_cardinality
        kins = sum(s(j,:));
        %ksum = sum(adj_mat(s_index(j),:));
        %kouts = ksum - kins;
        kins_sum = kins_sum + kins;
        %kouts_sum = kouts_sum + kouts;
        ec=kins_sum;
    end
    cf_s = (ec*1.0+2*s_cardinality)/(s_cardinality);
%     cf_s = (ec*1.0)/(s_cardinality);
    cf = cf + cf_s;
%     cf_s = (ec*1.0)/(2*ll);
%     cf = cf + cf_s;
end
%cf =0.5-cf;
cf=2*length(adj_mat)-cf;
end

