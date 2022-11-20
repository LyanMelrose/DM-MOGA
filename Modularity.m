function Q = Modularity(solution,AdjacentMatrix)
%% 获取无向无权网络的某种划分的模块度
%% 模块度的提出和计算方法见如下文献:
%% 简化版求Q
degree=single(sum(AdjacentMatrix,2));
edge_num=sum(degree)/2;
solution=single(solution);
m=max(solution);
Q=0;
for i=1:m
    if length(find(solution==i))>0
    Community_AdjMatrix=logical(AdjacentMatrix(find(solution==i),find(solution==i)));
    Degree_Matrix=single(degree(find(solution==i)))*single(degree(find(solution==i))')/(2*edge_num);
    Q=Q+sum(sum(Community_AdjMatrix-Degree_Matrix,1));
    end
end
clear Community_AdjMatrix Degree_Matrix;
Q=Q/(2*edge_num);



end
