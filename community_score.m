% function [ cs ] = community_score(adj_mat,clu_assignment,ll)
% % Compute the community score for each partition.
% % adj_mat: the adjacency matrix of the network.
% % clu_assignment: the cluster label vector.
% 
% r = 1;  % a parameter in the formula of community score.
% cs = 0;
% clu_num = max(clu_assignment);
% de=0;
% for i = 1:clu_num
%     s_index = find(clu_assignment == i);
%     s = adj_mat(s_index,s_index);
%     s_cardinality = length(s_index);
%     kins_sum = 0;
%     kouts_sum = 0;
%     for j = 1:s_cardinality
%         kins = sum(s(j,:));
%         ksum = sum(adj_mat(s_index(j),:));
%         kouts = ksum - kins;
%         kins_sum = kins_sum + kins;
%         kouts_sum = kouts_sum + kouts;
%         de=kouts_sum;
%     end
%     cf_s = de*1.0/(s_cardinality);
%     cs = cs + cf_s;
% end
% cs =cs;
% end

function [sum_v]=community_score(maxtrix,clu_assignment,ll)
global Node1 degree
a=1;
ec=0;
t=length(maxtrix);
for m=1:t
    B(m)=Node1(m).e(1);
    
end
[~,C]=sort(B);
sum_v=0;
clu_num=max(clu_assignment);
 for i=1:clu_num
    vertice=0;
    s_index=find(clu_assignment==i);
    sum1=0;
    sum2=0;
    A=[];
    for j=1:length(s_index)
        ec_1=0;
        for k=j+1:length(s_index)
            ec_1=ec_1+maxtrix(s_index(j),s_index(k));
        end
        ec=maxtrix(s_index(j),s_index(j));
        sum1=sum1+ec+ec_1;
        cc=find(C==s_index(j));
        A=[A Node1(cc).e];
        
        vertice=size(A,2);
       
    end
    sum2=sum(degree(A),2)/2;
        sum1=sum2-sum1;
    sum_v=sum1*2/vertice+sum_v;
        
end
sum_v=sum_v;
end
