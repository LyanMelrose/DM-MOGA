function Q = Qodularity(solution,AdjacentMatrix)
%% ��ȡ������Ȩ�����ĳ�ֻ��ֵ�ģ���
%% ģ��ȵ�����ͼ��㷽������������:
%% �򻯰���Q
degree=single(sum(AdjacentMatrix,2));
edge_num=sum(degree)/2;
solution=single(solution);
m=max(solution);
Q=0;
for i=1:m
    if length(find(solution==i))>0
        Community_AdjMatrix=logical(AdjacentMatrix(solution==i,solution==i));
        Degree_Matrix=single(degree(solution==i))*single(degree(solution==i)')/(2*edge_num);
        Q=Q+sum(sum(Community_AdjMatrix-Degree_Matrix,1));
    end
end
clear Community_AdjMatrix Degree_Matrix;
Q=Q/(2*edge_num);



end
