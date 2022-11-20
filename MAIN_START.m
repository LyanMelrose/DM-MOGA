function MAIN_START(path,real_path,name)
starttime = clock;

tic;

global   numVar  max_gen  M  pc pg pm
global   strNetwork BinaryAdj
global Thirty_Run_maxWW  AdjMatrix t3   C_result
global original_Matrix    label    gene_id_symbol    GT_gene    sim_Matrix


strNetwork=name;
pg=1;
max_gen=100;
M = 2;
pm=0.1;
pc=1;
runtimes=1;
All_Run_maxWW=[];   %存储三十次实验，每次WW值最大的结果的Q值和NMI值
All_Run_maxWW_Cres=[];   %存储三十次实验，每次WW值最大的结果的Q值和NMI值

%% 读入邻接矩阵或边表
for mg=1:runtimes
    load(path);
    load('MalaCards_all_connected_genes.mat');
    AdjMatrix=after_PPI_Network;
    BinaryAdj=AdjMatrix~=0;
    degree=uint16(sum(BinaryAdj,1));
%     % 从邻接矩阵中去掉所有孤立的点
%     index=find(degree==0);
%     AdjMatrix(index,:)=[];
%     AdjMatrix(:,index)=[];
    numVar=length(AdjMatrix);
    
    %% 获取合并后的团及相应的信息
    t3=clock;
    % 预缩减网络
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/6/16%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     [Clique,Matrix,degree1,edgeslist]=local_expansion(AdjMatrix);
    Clique=local_expansion(AdjMatrix);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2020/12/27%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MCODE vertex weighting
    disp('==========VertexWeighting==========')
    Weights=VertexWeighting(BinaryAdj);
    % 得到标号表示的社团划分结果
    Pre_Community = Clique2CIndex(Clique);

    Branch_Set=ones(numVar,1);
    % 初始化修剪枝杈后的社团划分
    New_Community=Pre_Community;
    for i=1:length(Clique)
        i
        % 对于第i个社团
        C_node=cell2mat(Clique(i));
        % 如果社团中大于等于3个节点，在New_Community和Branch_Set中将这些节点的对应位置标记为0
        if length(C_node)>=3
%             Branch_Set(C_node,1)=0;
%             New_Community(C_node,1)=0;
%         else
            % 提取第i个社团的节点及其权重，对该社团进行修剪,将剪掉的节点在New_Community和Branch_Set中标记为0
            [New_Clique,Branch_node]=Cut_branch(AdjMatrix,C_node,Weights(C_node),i,length(Clique));
            Branch_Set(Branch_node,1)=0;
            New_Community(Branch_node)=0;
            % 提取改变社团的节点，以及被修剪的节点
            ChangeID=New_Clique(New_Clique(:,2)~=New_Clique(:,3),1);
            New_Community(ChangeID,1)=New_Clique(New_Clique(:,2)~=New_Clique(:,3),3);
        end
    end
    Cut_Clique=CIndex2Clique(New_Community);
    Cut_Clique(cellfun(@isempty,Cut_Clique))=[];
    % 找到所有被修剪的节点，将单个节点分别作为一个局部社团加入Clique列表中
    Branch_node=find(Branch_Set==0);
    Clique_count=length(Cut_Clique);
    for i=1:length(Branch_node)
        Cut_Clique{Clique_count+1}=Branch_node(i);
        Clique_count=Clique_count+1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/6/16%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 矩阵压缩时怎么处理边的权重很重要。
    Matrix=flod(BinaryAdj,AdjMatrix,Cut_Clique);
    Matrix=single(Matrix);
    
    % edgeslist中包含大量空值行，可能的解决方法有二：
    % 一、矩阵压缩的时候错误，修改bug；
    % 二、这是正常现象，修改后续代码以适应；
    % 因此关键在判断这是否为正常现象。确认孤立的局部社团中的节点是否真的没有外连边。
    edgeslist = edges_list(Matrix,length(Cut_Clique),Cut_Clique);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/1/4%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    edgeslist([edgeslist(:).n]==0)=[];

    % 指定'UniformOutput'的值为false，cellfun()返回元胞数组类型的结果。
    degree1=cell2mat(cellfun(@(S)sum(degree(1,S)),Cut_Clique,'UniformOutput',false));


   %% 文件数据保存及进化得到parate面
    root = sprintf('results1/%s/ParetoFront',strNetwork);  %%建个文件保存实验数据
    if ~isdir(root) %判断路径是否存在
        mkdir(root);
    end
    root = sprintf('results1/%s/metrics',strNetwork);
    if ~isdir(root) %判断路径是否存在
        mkdir(root);
    end
    
    %这里的Thirty_Run其实存储的是每次实验的结果，而不是三十次实验的结果都有，三十次饰实验的总结果见All_Run。
    Thirty_Run_maxWW=0;   %存储本次实验Q值最大的结果的Q值，NMI值和RI值
    
    fprintf('第%s次运行\n',num2str(mg));
    
    [C_result,]=TMOEAD(Matrix,Cut_Clique,BinaryAdj,mg,degree1,starttime,edgeslist);
    
    %将本次实验的Thirty_Run汇总在All_Run中。
    All_Run_maxWW=[All_Run_maxWW;Thirty_Run_maxWW];
    All_Run_maxWW_Cres=[All_Run_maxWW_Cres;C_result];
    
    fprintf('\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/8/16%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,index1]=max(All_Run_maxWW(:,1));    %找到最大的 Weighted Wscore 索引

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/8/16%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 选择基因个数最大的10个模块
Cluster_F=tabulate(All_Run_maxWW_Cres(index1,:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/8/17%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Selected_Modules=sortrows(Cluster_F,-2);
Selected_Modules(Selected_Modules(:,1)==0,:)=[];
Selected_Modules=Selected_Modules(1:10,1);
% Selected_Modules=Cluster_F(Cluster_F(:,2)>100,1);
%Selected_Modules(Selected_Modules==0)=[];

% 依次得到每个被选中的模块的基因并输出，第一列为 gene symbol，第二列为模块标号（将原始标号改为1234顺排）
fid=fopen(['results1/' name '/metrics/' strNetwork '_gene.txt'],'wt');
for Cid=1:length(Selected_Modules)
    Ci_gene=gene_id_symbol(All_Run_maxWW_Cres(index1,:)==Selected_Modules(Cid));
    fprintf(fid,'%s\n',['Cluster Index: ' num2str(Cid)]);
    out=Ci_gene;
    [m,n]=size(out);
    for i=1:m
        for j=1:n
            x=out{i,j};
            fprintf(fid,'%s\t%s\n',x,num2str(Cid));
        end
    end
end
fclose(fid);

fprintf('Thirtyrun_results:\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/8/16%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('maxWW =    %g\n',All_Run_maxWW(index1,1));
path = sprintf('results1/%s/metrics/MyRMOEA_%s_minWW_%d.txt',strNetwork,strNetwork);
savedata1(path,[Selected_Modules' 0 0 0 All_Run_maxWW_Cres(index1,:)]);

clear;
