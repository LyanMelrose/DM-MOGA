
function [clique] =find_overlap_node(clique,over_node,AdjMatrix)
h=[];
over_node=unique(over_node);
% 对于每个重叠节点
for i=1:length(over_node)
    % cellfun(@function(),)可以利用一些定义的函数批量处理cell2mat产生的子矩阵，每个元胞为一个子矩阵。
    % 得到每个局部社团中第i个重叠节点的个数
    A=cell2mat(cellfun(@(s)length(find(s==over_node(i))),clique,'UniformOutput',false));
    B=find(A>0);
    %h=[h B];
    % a为重叠节点i的标号
    a=over_node(i);
    edges=[];
    % 对于每一个包含重叠节点i的局部社团，，，
    for j=1:length(B)
        % 统计局部社团中节点的内部连接度
        degree_C=sum(AdjMatrix(clique{B(j)},clique{B(j)}));
        index=find(clique{B(j)}==a);
        % 统计局部社团中与重叠节点i相连的边
        edges=[edges degree_C(index)];
    end
    %% 仅保留连边数最多的局部社团中的重叠节点i，其余的局部社团中全部删掉这个节点
    [~,index]= max(edges);
    % 从B中删除连边数最多的那个局部社团，并升序排序
    B=setdiff(B,B(index));
    % 删掉其他局部社团中的重叠节点i
    for k=1:length(B)
        clique{B(k)}=setdiff(clique{B(k)},a);
    end
    % end of over_node(i)
end
a=[];
% 重新组织存储局部社团
for i=1:length(clique)
    a=[a clique{i}];
end
length(a)

% community=clique(unique(h));
% clique(unique(h))=[];
% 
% 
% 
% 
% t=length(community);
% for i=1:t
%     
%    for j=1:t
%        if i~=j
%         A=community{i};
%         B=community{j};
%         t1=length(B);
%        
%         t2=length(A);
%         if t2<t1
%             t1=t2;
%         end
%         M=intersect(A,B);
%         t3=length(M);
%         if t3/t1>=0.5
%             community{i}=unique([community{i} community{j}]);
%         end
%        end
%    end
% 
% end
%   C=[];
% for i=1:t
%   A=community{i};
%    for j=1:t
%        if i~=j
%         B=community{j};
%         M=intersect(A,B);
%         t1=length(B);
%         t2=length(A);
%         if t2<t1
%             t1=t2;
%         end
%         
%         t3=length(M);
%         if t3>0
%             stop=1;
%         end
%         
%         if t3/t1>=0.5
%             community{i}=unique([community{i} community{j}]);
%             community{j}=[];
%         else
%             C=[C M];
%         end
%         end
%        
%    end
% end
% C=unique(C);
% m=1;
% for i=1:t
%     if ~isempty(community{i})
%         Clique{m}=setdiff(community{i},C);
%         m=m+1;
%     end
% end
% c=length(C);
% if c>0
% for i=1:c
%     Clique{m}=C(i);
%     m=m+1;
% end
% k=length(clique);
% for i=1:length(Clique)
%     clique{k+1}=Clique{i};
%     k=k+1;
% end
% end
% 
% 
% 
% 
% 
% % t=length(clique);
% % Q=[];
% % for i=1:t
% %  Q=sort([Q clique{i}]);
% % end
% end
%     
%     
%     
%     
% 
% 
% 
% 
