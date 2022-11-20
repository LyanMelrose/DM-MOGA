function [ cs ] = community_score(adj_mat,clu_assignment,ll)
% Compute the community score for each partition.
% adj_mat: the adjacency matrix of the network.
% clu_assignment: the cluster label vector.

r = 1;  % a parameter in the formula of community score.
cs = 0;
clu_num = max(clu_assignment);
for i = 1:clu_num
    s_index = find(clu_assignment == i);
    s = adj_mat(s_index,s_index);
    s_cardinality = length(s_index);
    vs = sum(sum(s));
    ms_temp = 0;
    for j = 1:s_cardinality
        kins = sum(s(j,:));
        u = kins/s_cardinality;
        ms_temp = ms_temp + u^r;
    end
    ms = ms_temp/s_cardinality;
    score_s = ms * vs;
    cs = cs + score_s;
end 
cs=-cs;
end

