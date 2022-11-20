function [ Clique,f,idealp,degree1,edgeslist,V,Matrix]=change_evolution(chromosomes,V,CLique,BinaryAdj,degree1,edgeslist,matrix,ori_matrix)
% function [ Clique,f, degree1, edgeslist,V,Matrix]=change_evolution(chromosomes,V,CLique,BinaryAdj,ll,degree1,edgeslist,matrix)
global numVar
% g 为目前网络中节点的个数（包括缩减后的局部社团）
g=V;      %用来控制最后无法合并的情况
% 群体中染色体的个数
popsize=size(chromosomes,1);
% 目前网络中节点的个数（包括缩减后的局部社团）

% 取群体中所有染色体的两个适应度函数值
BB=chromosomes(:,[V+1 V+2]);
% 对所有染色体进行非支配排序
B=P_sort(BB,'all');
% 得到 pareto 第一前沿的成员的索引，将它们作为精英个体
B=find(B==1);
ParetoFront=chromosomes(:,1:V);% 提取所有染色体的内容，不包括适应度函数值。
t1=clock;
%%%对个体团内的点进行震荡，找出团内变化的点能对Q产生正面作用。
% 对于Pareto前沿中的所有个体
for i=1:size(ParetoFront,1)
    % 将其解码为各个模块标号
    ParetoFront1(i,:) = decode1(ParetoFront(i,1:V),CLique);
    C_num(i,1)=max(ParetoFront1(i,:));   %%统计所有个体的模块划分结果中模块的个数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/8/31%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    WW(i,1) = Weighted_Wscore(ParetoFront1(i,:),ori_matrix);%%计算每个个体的加权 W 值
end


% 找到精英个体中划分错误的节点
% 返回修改后的精英个体模块划分结果 ParetoFront1(B,:),以及错误节点的索引 remove。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/8/16%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ParetoFront1(B,:),remove]=find_error(ParetoFront1(B,:),ori_matrix,CLique,WW);
% [ParetoFront1(B,:),remove]=find_error(ParetoFront1(B,:),matrix,CLique,WW);
%%%标签局部扩散----------------
Clique=find_merge(WW,C_num,chromosomes,CLique);
% 生成新的局部社团。
Clique=find_clique(Clique,remove,numVar);
Clique=Clique';
% 得到所有局部社团更新后，局部社团的个数 V
V=length(Clique);
% 如果局部社团个数发生改变
if V~=g
    % 进行邻接矩阵网络缩减
    Matrix=find_Matrix(ori_matrix,Clique);
    % 更新边表信息和节点度信息
    edgeslist = edges_list(Matrix,V,Clique);
    degree=sum(BinaryAdj,1);
    degree1=cell2mat(cellfun(@(S)sum(degree(S)),Clique,'UniformOutput',false));
    
    %%%继承前面30代的模块划分信息，即虽然进行了网络缩减，生成了新的局部社团，但其中节点所属模块保持不变%%%%
    % 对于每个个体
    for i=1:popsize
        % 将前30代中个体 i 所获每个局部社团的模块索引赋给该局部社团的所有节点的新模块索引
        for j=1:length(Clique)
            ParetoFront1(i,Clique{j})=ParetoFront1(i,Clique{j}(1));
        end
        % inherit为继承函数，将前面修改后的 ParetoFront1 带入，但返回的向量是编码后的结果，即chromosomes
        f(i,1:V)=inherit(ParetoFront1(i,:),Clique);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/7/30%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        f(i,V + 1: V+2) = evaluate_objective(f(i,1:V),degree1,Matrix,Clique);
        %         f(i,V + 1: V+2) = evaluate_objective(f(i,1:V),ll,degree1,Matrix,Clique);
    end
    % 更新适应度函数全局最小值
    idealp = min(f(:,V+1:V+2));
else      % 如果局部社团个数未发生改变
    % 不更改chromosomes，CLique，以及matrix
    f=chromosomes;
    idealp = min(f(:,V+1:V+2));
    Clique=CLique;
    Matrix=matrix;
end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%函数部分%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ParetoFront1,remove]=find_error(ParetoFront1,AdjMatrix,CLique,W)
M=[];
% V 为目前网络中节点的个数（包括缩减后的局部社团）
V=length(CLique);
% numVar 为未缩减网络中的节点个数
numVar=size(AdjMatrix,1);

% 对于每个个体
for in=1:size(ParetoFront1,1)
    % 初始化 change_node
    change_node{in}=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/8/16%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % label，提取第 in 个个体的模块划分结果
    label= ParetoFront1(in,1:numVar);
%     label= ParetoFront1(in,1:numVar);
    
    %% %% 后处理：2种方式-最大的Weighted Wscore label = ParetoFront1(in,1:numVar);
    % 对于每个节点（包含局部社团）
    for i=1:V
        % 如果局部社团中包含大于三个节点
        if(length(CLique{i})>=3)
            % 对于局部社团中的每个节点
            for j=1:length(CLique{i})
                % 局部社团中保存了所含节点的索引，提取第 j 个节点的索引 m
                m=CLique{i}(j);
                % 提取第 j 个节点所属模块的编号 k
                k=label(m);
                % index_last，找到模块 k 中所有节点的索引
                index_last=find(label==k);
                % 检查节点 m 是否包含于旧的局部社团中
                A=setdiff(index_last,m);
                % neighbors，提取节点 m 的邻居
                neighbors=find(AdjMatrix(m,:)~=0);
                % Max()为自定义函数，用于寻找包含邻居个数最多的模块的id
                com_max=Max(neighbors,label);
                % com_max 返回节点 m 的邻居所属模块，该模块包含其邻居的个数最多。
                % 以下语句使得节点 m 也属于这一模块。
                ParetoFront1(in,m)=com_max;
                % 如果节点 m 的所属模块发生了改变
                if com_max~=k
                    % 将 m 标记为第 in 个局部社团的 change node，并保存。
                    change_node{in}=[change_node{in} m];
                end
            end
        end
    end
    % ParetoFront1 已经在第 101 行进行了更新，计算新的 Weighted Wscore
    WW(in) = Weighted_Wscore(ParetoFront1(in,1:numVar),AdjMatrix);
    % 如果change node中的节点修改后能够获得更大的 Weighted Wscore 值，提取这个个体
    if WW(in)<W(in)                                                       %有疑问
        M=[M in];
    end
end
error_node=[];
% 将M中的个体的change node 返回，作为划分错误的节点，以及函数返回值。
for i=1:length(M)
    error_node=[error_node change_node{M(i)}];
end
remove=unique(error_node);

end


function Clique=find_merge(W,C_num,chromosomes,CLique)
%需要两个指标越大越好
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/8/16%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f(:,1)=W;
% f(:,1)=-Q;

f(:,2)=-C_num;
V=length(CLique);
FrontValue = P_sort(f,'all');   %front层
B=find(FrontValue==1);
K=chromosomes(B,1:V);
[a,~]=size(K);
MM=[];
for i=1:a
    A=decode(K(i,:));
    MM=[MM;A];
end
visited=zeros(1,V);
times=1;
for i=1:V
    if visited(i)==0
        distance_i=zeros(1,V);
        for j=i+1:V
            if visited(j)==0
                % 统计局部社团i与j不同的位点个数，a应该是个体的数量
                A=MM(:,i)-MM(:,j);
                index_0=find(A==0);
                distance_i(j)=length(index_0)/a;
            end
        end
        % 找到在大部分个体中相同的局部社团
        F=find(distance_i>0.9);
        visited(F)=1;
        % times有什么作用呢
        Clique{times}=cell2mat(CLique([F i]));
        if length(Clique{times})>0
            times=times+1;
        end
    end
end
end



function Clique=find_clique(Clique,erase_node,numVar)
t=length(Clique);
for i=1:t
    Clique{i}=setdiff(Clique{i},erase_node);
end
R=[];
for i=1:t
    R=[R Clique{i}'];
end
R=setdiff(1:numVar,R);
for i=1:length(R)
    t=t+1;
    Clique{t}=R(i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/8/15%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Clique=Clique(find(cell2mat(cellfun(@(S)length(S),Clique,'UniformOutput',false))~=0));  %%剔除长度为0的子团！

[~,rank] = sort(cell2mat(cellfun(@(s)s(1),Clique,'UniformOutput',false)));                %对元胞进行排序
Clique= Clique(rank);
end

function Matrix=find_Matrix(Matrix,Clique)
numVar=size(Matrix,1);
% degree=sum(Matrix,1);
D=1:numVar;
for i=1:length(Clique)
    C=Clique{i};
    vertex_min=find(D==C(1));
    while length(C)>1
        j=2;
        vertex_max=find(D==C(2));
        if length(vertex_max)==0 %如果在2758个压缩后节点中找不到第二个节点，那么不可以压缩
            CD_node=find(Clique{i}==C(j));
            Clique{i}(CD_node)=[];%将第二个节点从局部社团中除去
            C(j)=[];
        else %如果在2758个压缩后节点中能够找到第二个节点，那么不可以压缩
            Matrix(:,vertex_min)=Matrix(:,vertex_min)+Matrix(:,vertex_max);
            Matrix(vertex_min,:)=Matrix(vertex_min,:)+Matrix(vertex_max,:);
            Matrix(vertex_min,vertex_min)=Matrix(vertex_min,vertex_min)-Matrix(vertex_min,vertex_max)-Matrix(vertex_max,vertex_max);
            Matrix(vertex_max,:)=[];
            Matrix(:,vertex_max)=[];
            index=find(D==C(2));
            D(index)=[];
            C(j)=[];
        end
    end
end
end

function f=inherit(ParetoFront1,Clique)

RANK=cell2mat(cellfun(@(S)S(randi(length(S))),Clique,'UniformOutput',false));
% 得到每个局部社团的模块索引
RANK=ParetoFront1(RANK);
% 考虑到有些局部社团属于同一个模块，因此下面的循环依次选出每个包含局部社团的模块，
for j=1:max(RANK)
    % 提取属于第 j 个模块的局部社团的索引。
    O=find(RANK==j);
    % 如果 O 非空
    if length(O)>0
        % 将 O 的第一个元素转而放在 O 的最后，
        O=[O O(1)];
        O(1)=[];
        % 因为 f 为邻居表编码后的结果，因此下列语句使得属于同一个模块的局部社团互连为环。
        f(find(RANK==j))=O;
    end
end

end