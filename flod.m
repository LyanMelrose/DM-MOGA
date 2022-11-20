function AdjMatrix=flod(BinaryAdj,matrix,clique)

% 得到所有节点的度
degree=sum(BinaryAdj,1);
% 得到局部社团的个数
len=length(clique);
% 得到节点个数
numVar=length(BinaryAdj);
% 初始化压缩后的邻接矩阵
AdjMatrix=zeros(len,len);
% 初始化节点社团标号向量，
label=zeros(1,numVar);
for i=1:len
    label(clique{i})=i;
end

% 节点的邻居的标签
for i=1:numVar
    D{i}=label(find(BinaryAdj(i,:)==1));
end

% 对于每个局部社团
for i=1:len
    % 提取局部社团i中的节点
    a=clique{i};
    % C为局部社团i的邻居局部社团
    C=[];
    for j=1:length(a)
        C=[C D{a(j)}];
    end
    C=unique(C);
    neighbor=unique(C);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/8/16%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 团i内部的边数
    AdjMatrix(i,i)=sum(sum(matrix(clique{i},clique{i})))/2;
    %     AdjMatrix(i,i)=sum(sum(BinaryAdj(clique{i},clique{i})))/2;
    
    % 对于局部社团i的每个邻居社团a，计算二者之间的连接边数
    for j=1:length(C)
        a=C(j);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2020/12/23%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        AdjMatrix(a,i)=sum(sum(matrix(clique{i},clique{a})));
%         AdjMatrix(a,i)=sum(sum(BinaryAdj(clique{i},clique{a})));
        AdjMatrix(i,a)=AdjMatrix(a,i);
    end
    AdjMatrix(i,i)=AdjMatrix(i,i)/2;
end

end