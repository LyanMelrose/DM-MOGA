function C_Result=TMOEAD(Matrix,Clique,BinaryAdj,mg,degree1,starttime,edgeslist)
%% 主代码部分
% 网络缩减操作仅与社团拓扑相关，各节点所属模块的搜索是遗传算法的主要目的，二者代码基本是分开的。

global popsize numVar  max_gen  M niche V pc pm
global strNetwork weights neighbors
global Thirty_Run_maxWW   Adj_mat   t3
global original_Matrix    label    gene_id_symbol    GT_gene    AdjMatrix

%% 进化前期处理--kk表示目标函数个数--ideal表示理想点--weight，neighbours表示权重和邻居--Chromosomes表示种群
kk=2;
popsize=50;
niche=40;


V=size(Matrix,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/8/2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idealp = Inf*ones(1,M);
% idealp = -Inf*ones(1,M);

% 初始化种群中个体的权重和邻居（应该是多目标分解的一步）
[weights neighbors] = init_weight(popsize, niche);
popsize = size(weights,1);
index=[];

% 初始化种群
chromosomes = initialize_variables(popsize,edgeslist,M,degree1,Matrix,Clique);
% chromosomes = initialize_variables(popsize,edgeslist,M,ll,degree1,Matrix,Clique);



%% 迭代进化
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/7/30%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stop=1;

Save_chromosome_sub_FF=zeros(max_gen,2);
% 将每一代中所选模块的各个评价指标值依次分别存储为三个矩阵
% 如果所选模块个数大于10，则增加列数
% 设初始选择模块的个数为0
CNUM=0;
Module_ACC=[];
Module_AUC=[];
Module_P=[];
All_Upper_Idx=zeros(max_gen,7);
for Gene = 1 : max_gen
    t1=clock;
    % max_gen=100，因此仅当 Gene=30 或 60 时执行进化过程中的网络缩减过程
%     if mod(Gene,30)==0&&Gene<=60
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/7/30%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %         stop=stop+1;
%         
%         % [ Node3,chromosomes, degree1, edgeslist,V,fault_node]= Fault_tolerance(Node3,chromosomes, degree1, edgeslist,V,fault_node,AdjMatrix);            %%加相关处理，处理错误的点，更新adj_mat与edgeslist，degree1，V
%         Draw(chromosomes,V);
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/7/30%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         [ CLique,chromosomes,idealp,degree1,edgeslist,V,Matrix] = change_evolution(chromosomes,V,Clique,BinaryAdj,degree1,edgeslist,Matrix,AdjMatrix);
%         %         [ Clique,chromosomes, degree1, edgeslist,V,Matrix] = change_evolution(chromosomes,V,Clique,BinaryAdj,ll,degree1,edgeslist,Matrix);
%         Clique=CLique;
%     end
    
    % 生成下一代个体。
    % 对于当前种群中的每个个体
    for i=1:popsize
        % 提取第 i 个个体的邻居
        i_neighbor_index = neighbors(i,:);
        % 提取其邻居的权重
        i_neighbor_weight = weights(i_neighbor_index,:);
        % 提取其邻居的染色体
        i_neighbor_chromosome = chromosomes(i_neighbor_index,:);
        % 提取其邻居的适应度函数值
        old_objective = i_neighbor_chromosome(:,V+1:V+M);
        
        % 遗传进化算子的操作函数
        % 这里修改了染色体，却没有修改clique,是否应该与change_evolution()的返回值相同？
        % 即，[ Clique,f,idealp,degree1,edgeslist,V,Matrix]
        % 不太对，虽然进行了模块划分的转变，但未进行网络缩减。
        %% 通过交叉、变异遗传算法算子获得两个新的染色体
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/8/2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        new_chromosomes = genetic_op(i_neighbor_chromosome,chromosomes,edgeslist,pc,pm,V,M,degree1,Matrix,Clique);
        %         new_chromosomes = genetic_op(i_neighbor_chromosome,chromosomes,edgeslist,pc,pm,V,M,ll,degree1,Matrix,Clique);
        
        % 更新全局最小适应度函数值，update idealp
        idealp = min(min(new_chromosomes(:,V+1:V+M)),idealp);


        % update the neighbours
        % 选择对于当前子问题最好的解
        chromosome_sub = sort_sub(weights(i,:),new_chromosomes,idealp,V,M);
        Save_chromosome_sub_FF(Gene,:)=chromosome_sub(V+1:V+M);
        
        %-------------------
        new_subobj = subobjective_te(i_neighbor_weight, chromosome_sub(V+1:V+M), idealp);
        old_subobj = subobjective_te(i_neighbor_weight, old_objective, idealp);
        jj=0;
        iposition = find(new_subobj<old_subobj);
        jj=length(iposition);
        if jj>kk         %%挑取2个进行替代
            kk1=floor(jj*rand(kk,1))+1;%%rand（2，1）产生随机的2行一列随机数
            iposition1=iposition(kk1);
        else
            iposition1=iposition;
        end
        chromosomes(i_neighbor_index(iposition1)',:) = chromosome_sub(ones(length(iposition1),1),:);
        Problem='聚类问题';
        M=2;
        clc;
        fprintf('NSGA-II,第%2s轮,%5s,%2s维,已完成%4s%%,耗时%5s秒\n',num2str(1),Problem,num2str(M),num2str(roundn(Gene/max_gen*100,-1)),num2str(roundn(toc,-2)));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/8/10%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 将本次循环的最好个体解码为非压缩社团划分结果。
    D_chromosome_sub=decode1(chromosome_sub(1,1:V),Clique);
    % 在这里评价适应度函数值和p值和分类测度的关系
    C_Frequency=tabulate(D_chromosome_sub);
    % 选择模块大小最大的前十个模块
    Selected_Modules=sortrows(C_Frequency,-2);
    Selected_Modules(Selected_Modules(:,1)==0,:)=[];
    Selected_Modules=Selected_Modules(1:10,1);
    % 更新所选模块的个数
    CNUM=10;
    % 对于每个被选中的模块，分别评价并记录指标值
    ACC=zeros(size(Selected_Modules'));
    AUC=zeros(size(Selected_Modules'));
    P=zeros(size(Selected_Modules'));
    for SC_idx=1:CNUM
        [Indicators,Avg_AUC] = FiveFold_SVM(original_Matrix(D_chromosome_sub==Selected_Modules(SC_idx,1),:)',label);
        ACC(SC_idx)=Indicators.accuracy;
        AUC(SC_idx)=Avg_AUC;
        
        gene_id_symbol_1=cellstr(gene_id_symbol);
        GT_gene_1=cellstr(GT_gene);
        K=length(intersect(gene_id_symbol_1,GT_gene_1));
        N=C_Frequency(Selected_Modules(SC_idx,1),2);
        P_value=hygepdf(1,length(gene_id_symbol),K,N);
        P(SC_idx)=P_value;
    end
    Gene
    FLE_It = FLE(D_chromosome_sub',BinaryAdj);
    All_Upper_Idx(Gene,:)=[Gene mean(ACC) mean(AUC) mean(P) mean(FLE_It(Selected_Modules)) idealp];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/8/15%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 画图,(1,1)为ACC，(1,2)为F1，(2,1)为AUC，(2,2)为PP
% figure(1)
% cmap=colormap();
% Col=64/(CNUM+1);
% plot((1:length(Module_ACC))',Save_chromosome_sub_FF(Module_ACC(:,1)~=0,1),'o','MarkerSize',4,...
%     'MarkerFaceColor',cmap(int8(1*Col),:),'MarkerEdgeColor',cmap(int8(1*Col),:));
% plot((1:length(Module_ACC))',Save_chromosome_sub_FF(Module_ACC(:,1)~=0,2),'o','MarkerSize',4,...
%     'MarkerFaceColor',cmap(int8(2*Col),:),'MarkerEdgeColor',cmap(int8(2*Col),:));
% plot((1:length(Module_ACC))',Module_ACC(Module_ACC(:,1)~=0,1),'o','MarkerSize',4,...
%     'MarkerFaceColor',cmap(int8(3*Col),:),'MarkerEdgeColor',cmap(int8(3*Col),:));
% plot((1:length(Module_ACC))',Module_AUC(Module_ACC(:,1)~=0,1),'o','MarkerSize',4,...
%     'MarkerFaceColor',cmap(int8(4*Col),:),'MarkerEdgeColor',cmap(int8(3*Col),:));
% plot((1:length(Module_ACC))',Module_ACC(Module_P(:,1)~=0,1),'o','MarkerSize',4,...
%     'MarkerFaceColor',cmap(int8(5*Col),:),'MarkerEdgeColor',cmap(int8(3*Col),:));
% hold on;
% xlabel('Iteration');
% ylabel('Evaluation Criterion');
% title('2D Nonclassical multidimensional scaling','FontSize',10);
% hold off;


%% 进化收敛保存parate面上的
times=etime(clock,t3);
Population=chromosomes;
%整理得到 pareto 前沿
ParetoFront=unique(Population,'rows');
path = sprintf('results1/%s/ParetoFront/MODPSO1_%s_PF%d.txt',strNetwork,strNetwork,mg);  %%种群进化最前沿数据
%% save_metrics
%ParetoFront1即为社团划分结果，其中每行代表一个Pareto最优解，每列代表一个节点所属的社团。
WW=zeros(size(ParetoFront,1),1);
for i=1:size(ParetoFront,1)
    ParetoFront1(i,1:numVar) = decode1(ParetoFront(i,1:V),Clique);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/8/31%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    WW(i,1) = Weighted_Wscore(ParetoFront1(i,:),AdjMatrix);%%计算每个个体的加权 W 值
end

% metrics=ParetoFront;
[~,index1]=max(WW);
fprintf('Weighted-W值最大时候的社团数目：%d\n',length(unique(ParetoFront1(index1,:))));
C_Result=ParetoFront1(index1,:);

disp(strNetwork);
fprintf('第%d次实验结果\n',mg);

savedata1(path,ParetoFront);
Thirty_Run_maxWW=WW(index1,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/8/25%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path = sprintf('results1/%s/metrics/%s_All_Upper_Idx_%d.txt',strNetwork,strNetwork,mg);
% 导出 All_Upper_Idx
savedata1(path,All_Upper_Idx);
% path = sprintf('results1/%s/metrics/%s_PF_FitnessValue%d.txt',strNetwork,strNetwork,mg);
% savedata1(path,metrics(:,V+1:V+M));


%应当是实现的最后一个修正策略，执行完成后上面三条语句的Thirty_Run_maxWW和path中保存的数据会被替换（如果修正后的WW更大）
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/8/16%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
last_operate(strNetwork,BinaryAdj,mg,degree1,Clique);
% last_operate(A,strNetwork,BinaryAdj,mg,degree1,Clique);

total_time = etime(clock, starttime);
path_time = sprintf('results1/%s/metrics/total_time_%d.txt',strNetwork,mg);
savedata1(path_time,total_time);

fprintf('total time used %u\n', etime(clock, starttime));
end