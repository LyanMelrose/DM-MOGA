function clique=local_expansion(Matrix)

numVar=size(Matrix,1);
 % over_node返回局部社团间重叠的节点，clique返回所有找到的局部社团
[over_node,clique]=localexpansion(Matrix,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2020/6/4%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
over_node=uint16(over_node);
% find_overlap_node()并非做无用功
t1=clock;
[clique]=find_overlap_node(clique,over_node,Matrix);


% Q中统计所有被划分的节点
T=length(clique);
Q=[];
for i=1:length(clique)
    Q=[Q uint16(clique{i})];
end
% 将未被划入局部社团的点补在局部社团向量clique的最后
A=1:length(Matrix);
A=setdiff(A,Q);
for i=1:length(A)
    clique{T+1}=A(i);
    T=T+1;
end

clique(cellfun(@isempty,clique))=[];
t=length(clique);
% 对元胞进行排序，
% 为什么要排序？
[~,rank] = sort(cell2mat(cellfun(@(s)s(1),clique,'UniformOutput',false)));
clique= clique(rank);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2020/12/24%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% degree=sum(Matrix,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/1/2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 仅用来调试，验证clique确实有重叠划分的节点。
% count=zeros(numVar,1);
% for i=1:length(clique)
%     cc=cell2mat(clique(i));
%     count(cc,1)=count(cc,1)+1;
% end

t1=clock;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/6/16%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 以下省略大量非原创代码
% Matrix=flod(Matrix,clique);
% Matrix=single(Matrix);
% edgeslist = edges_list(Matrix,t,clique);
% degree1=cell2mat(cellfun(@(S)sum(degree(1,S)),clique,'UniformOutput',false));

end
