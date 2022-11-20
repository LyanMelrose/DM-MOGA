function f = evaluate_objective(x,degree1,matrix,clique)
global sim_Matrix

f=[];
% 解码
x1 = decode(x);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/9/3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x2=decode1(x,clique);
% f(1) = -Weighted_Wscore(x1,matrix);
% 将每个节点的行提取出来并求和，得到带权度，设带权度最大的是聚类中心
C_idx=0;
maxInnerW=[];
for i=1:max(x2)
    %找到属于簇i的下标值
    y=find(x2==i);
    [~,c_id]=max(sum(sim_Matrix(y,:),2));
    C_idx(i)=y(c_id);
end
for i=1:max(x1)
    %找到属于簇i的下标值
    y=find(x1==i);
%     maxInnerW=[maxInnerW alpha*sum(matrix(y,y),[1 2])+(1-alpha)*length(y)];
    maxInnerW=[maxInnerW sum(matrix(y,y),[1 2])];
    CC_value(i)=CCIndex( matrix(y,y) );
end
Matrix=1./sim_Matrix;
Matrix(isinf(Matrix))=0;

% f(1) = -sum(maxInnerW);
f(1) = -max(CC_value);
f(2) = DBIndex(Matrix,C_idx,x2);% Dunns(max(x2),sim_Matrix,x2);

% f = -Weighted_Wscore(x,matrix);
% f(2) = fun_c(matrix,x,degree1,clique);

end