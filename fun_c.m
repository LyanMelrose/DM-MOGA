% function [ cf ] = community_fitness_mod2( adj_mat,clu_assignment,ll )
% % Compute the community fittness for each partition.
% % adj_mat: the adjacency matrix of the network.
% % clu_assignment: the cluster label vector.
% 
% a = 1;  % a parameter in the formula of community fittness.
% cf = 0;
% ec=0;
% clu_num = max(clu_assignment);
% for i = 1:clu_num
%     s_index = find(clu_assignment == i);
%     s = adj_mat(s_index,s_index);
%     s_cardinality = length(s_index);
%     kins_sum = 0;
%     kouts_sum = 0;
%     for j = 1:s_cardinality
%         kins = sum(s(j,:));
%         %ksum = sum(adj_mat(s_index(j),:));
%         %kouts = ksum - kins;
%         kins_sum = kins_sum + kins;
%         %kouts_sum = kouts_sum + kouts;
%         ec=kins_sum;
%     end
%     cf_s = ec*1.0/(s_cardinality);
%     cf = cf + cf_s;
% end
% cf = -cf;
% end

function f=fun_c(Matrix,label,degree1,clique)

% a=1;
% ec=0;
% t=length(Maxtrix);
% 
% sum_v=0;
% cs_v=0;
% clu_num=max(lable);
%  for i=1:clu_num
%     vertice=0;
%     s_index=find(lable==i);
%     sum1=0;
%     for j=1:length(s_index)
%         ec_1=0;
%         for k=j+1:length(s_index)
%             ec_1=ec_1+Maxtrix(s_index(j),s_index(k));
%         end
%         ec=Maxtrix(s_index(j),s_index(j));
%         sum1=sum1+ec+ec_1;
% %        cc=find(C==s_index(j));
%      
%     end
%      nn=length(s_index);
%       B=[];
%       A=[];
%         for mm=1:nn
%             B=clique{s_index(mm)};
%         A=[A B];
%         
%         end
%         community_edge=sum(degree1(s_index));
%         edge=community_edge-sum1*2;
%         vertice=vertice+size(A,2);
%         A=[];
%     sum_v=sum1*2/vertice+sum_v;
%      cs_v=edge/vertice+cs_v;   
% end
% sum_v=-sum_v;
% cs_v=cs_v;
% f=[sum_v cs_v];

sum_v_in=0;
sum_v_out=0;
t=max(label);
C=0;
for i=1:t
    index=find(label==i);
    if length(index)>0
        m=Matrix(index,index);
        edges_in=(sum(sum(m))-sum(diag(m)))/2+sum(diag(m));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/8/31%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        edges_out=sum(sum(Matrix(index,label~=i)));
%         edges_out=sum(degree1(index))-edges_in*2;

        n=length(cell2mat(clique(index)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/8/31%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nc=length(cell2mat(clique))-n;
        sum_v_in=edges_in/n+sum_v_in;
        sum_v_out=edges_out/n+sum_v_out;
%         sum_v_in=edges_in/(n^2)+sum_v_in;
%         sum_v_out=edges_out/(n*nc)+sum_v_out;
%         sum_v_in=edges_in/n+sum_v_in;
%         sum_v_out=edges_out/n+sum_v_out;
        
        C=C+1;
    end
end

    f(1)=-sum_v_in*2;
    f(2)=sum_v_out;


end