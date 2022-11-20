function f= initialize_variables(N,el,M,degree1,matrix,Node3)
% function [f]= initialize_variables(N,el,M,ll,degree1,matrix,Node3)
% 这里的matrix是普通邻接矩阵缩减后的矩阵，未有
global idealp;
V=size(matrix,1);

K = M + V;

% N为个体个数，对于每个个体，保存在矩阵f中。
% 每行代表一个粒子，每列代表一个邻居局部模块的某个节点，元素中存储对应节点的邻居的标号
for i = 1 : N
    if rand<0.2
        f(i,1:V)=IGLP(V,el,Node3);
    else
        % el为edgeslist
        for j = 1 : V
            % 从第j个局部模块中随机选择一个节点的邻居标号
            f(i,j) = el(j).e(ceil(rand*el(j).n));
        end
    end
    %每行的最后两列存储目标函数值
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/8/31%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    f(i,V + 1: K) = evaluate_objective(f(i,1:V),degree1,matrix,Node3);
%     f(i,V + 1: K) = evaluate_objective(f(i,1:V),ll,degree1,matrix,Node3);

end

% 初始化为群中两个分别最小的目标函数值
idealp = min(f(:,V+1:K));
 
end
