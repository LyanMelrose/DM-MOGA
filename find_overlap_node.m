
function [clique] =find_overlap_node(clique,over_node,AdjMatrix)
h=[];
over_node=unique(over_node);
% ����ÿ���ص��ڵ�
for i=1:length(over_node)
    % cellfun(@function(),)��������һЩ����ĺ�����������cell2mat�������Ӿ���ÿ��Ԫ��Ϊһ���Ӿ���
    % �õ�ÿ���ֲ������е�i���ص��ڵ�ĸ���
    A=cell2mat(cellfun(@(s)length(find(s==over_node(i))),clique,'UniformOutput',false));
    B=find(A>0);
    %h=[h B];
    % aΪ�ص��ڵ�i�ı��
    a=over_node(i);
    edges=[];
    % ����ÿһ�������ص��ڵ�i�ľֲ����ţ�����
    for j=1:length(B)
        % ͳ�ƾֲ������нڵ���ڲ����Ӷ�
        degree_C=sum(AdjMatrix(clique{B(j)},clique{B(j)}));
        index=find(clique{B(j)}==a);
        % ͳ�ƾֲ����������ص��ڵ�i�����ı�
        edges=[edges degree_C(index)];
    end
    %% ���������������ľֲ������е��ص��ڵ�i������ľֲ�������ȫ��ɾ������ڵ�
    [~,index]= max(edges);
    % ��B��ɾ�������������Ǹ��ֲ����ţ�����������
    B=setdiff(B,B(index));
    % ɾ�������ֲ������е��ص��ڵ�i
    for k=1:length(B)
        clique{B(k)}=setdiff(clique{B(k)},a);
    end
    % end of over_node(i)
end
a=[];
% ������֯�洢�ֲ�����
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
