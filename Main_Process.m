%% 使用R包limma进行差异表达基因分析，同时计算 sim_Matrix，
% 导入 sim 相似性矩阵为 sim_Matrix，
% 导入 sim_genes_expression 的矩阵为 original_Matrix，样本导入为 sample_id，
% 基因导入为 sim_genes


%% 使用高斯核MI计算差异表达基因间的相互作用
% 导入 _ENid_BP_GSE18842_SoftFinal.mat 中的 original_Matrix，转至子文件夹./robince-gcmi-f14aae8/matlab
[gene_num,~]=size(original_Matrix);
Network=zeros(gene_num,gene_num);
for i=1:gene_num-1
    for j=(1+i):gene_num
        %I = kernelmi(Matrix(i,:),Matrix(j,:));
        %保证基因表达数据的行为样本，列为矩阵
        I = gcmi_cc(original_Matrix(i,:)',original_Matrix(j,:)');
        Network(i,j)=I;
        Network(j,i)=I;
    end
end

%min-max标准化
Net_2=(Network-min(Network(:)))/(max(Network(:))-min(Network(:)));
%Net_2为[0,1]
% Net_2=single(Net);
% Net_2(logical(eye(size(Net_2))))=1;
%Net_1为[-1，1]
% Net_1=(Net_2-0.5).*2;
% Net_1(logical(eye(size(Net))))=0;
%I = gcmi_cc(x, y)也可以用
Network=Net_2;

%% 通过 HPRD 的 PPI 筛选软阈值共表达网络的边
% 1. 运行 select_Diff_gene_p.m，将 sim_genes 的 ENTREZ ID 转化为 gene symbol,
% 2. 运行 select_HPRD_Edges_from_Adj.m，得到 PPI_Adj，
% 3. 用 PPI_Adj 来过滤上一步获得的 Network，
PPI_Adj=logical(PPI_Adj);
Net_2=zeros(size(Network));
Net_2(PPI_Adj)=Network(PPI_Adj);
Network=Net_2;

% 去掉 Network 中度为0的节点
non_zero_id=find(sum(Network)~=0);
Network=Network(non_zero_id,non_zero_id);
% 从 gene_id_symbol，gene_p，original_Matrix，sim_genes，sim_Matrix 中去掉度为 0 的节点。
% non_zero_id=find(sum(PPI_Adj)~=0);
gene_id_symbol=gene_id_symbol(non_zero_id,1);
gene_p=gene_p(non_zero_id,1);
sim_genes=sim_genes(non_zero_id,1);
original_Matrix=original_Matrix(non_zero_id,:);
sim_Matrix=sim_Matrix(non_zero_id,non_zero_id);
% 继而，Network 改名为 after_PPI_Network
% after_PPI_Network，gene_id_symbol，gene_p，original_Matrix，sim_genes，sim_Matrix，non_zero_id(node_id)
% 将被保存为 _ENid_BP_GSE19804_Final.mat


%% 转至Rstudio，使用WGCNA选取阈值
% 发现软阈值取为1，所以可以省掉这一步。


%% 执行My_RMOEA - 加权算法
% 导入 _ENid_BP_GSE19804_Final.mat，转至“My_RMOEA - 加权”文件夹
my_run


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/8/2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 实验一：通过五折交叉验证和SVM进行样本分类，以定量评价所选模块的分类效果。
% 分别取模块个数最大的十个取模块中的基因用来五折验证
% 导入 My_RMOEA-加权 获得的实验结果为 Cluster元胞数组（基因），同时导入选取的模块索引为 Cluster_id
for i=1:10
    Accuracy=[];
    Precision=[];
    Recall=[];
    F_measure=[];
    Average_AUC=[];
    % 对于每个社团计算 10 次
    for t=1:10
        Cluster_genes=Cluster(Cluster_id==i);
        selected_id=[];
        for j=1:length(Cluster_genes)
            bi_id=cell2mat(cellfun(@(x)Cluster_genes{j,1}==x,gene_id_symbol,'UniformOutput',false));
            selected_id=[selected_id;find(bi_id)];
        end
        [Indicators,Average_AUC(t,1)] = FiveFold_SVM(original_Matrix(selected_id,:)',label);
        Accuracy(t,1)=Indicators.accuracy;
        Precision(t,1)=Indicators.precision;
        Recall(t,1)=Indicators.recall;
        F_measure(t,1)=Indicators.F_measure;
        name='results1/_ENid_BP_GSE19188_PCC_Final_C/metrics/_ENid_BP_GSE19188_PCC_Final_C';
%         name='D:/_Documents/_Papers/Paper_3/Comparison_Methods/My_RMOEA/results1/_ENid_BP_GSE19188_Final_C/metrics/_ENid_BP_GSE19188_Final_C';
        write_Net=[name num2str(i) '_' 'Indicators.txt'];
        fid=fopen(write_Net,'w');%写入文件路径
        fprintf(fid,'%s\n',['ACC: ' num2str(mean(Accuracy))]);
        fprintf(fid,'%s\n',['Precision: ' num2str(mean(Precision))]);
        fprintf(fid,'%s\n',['Recall: ' num2str(mean(Recall))]);
        fprintf(fid,'%s\n',['F_measure: ' num2str(mean([Indicators.F_measure']))]);
        fprintf(fid,'%s\n',['G-mean: ' num2str(mean([Indicators.G_mean]))]);
        fprintf(fid,'%s\n',['Average_AUC: ' num2str(mean(Average_AUC))]);
    end
end
    
%% 实验二：通过MalaCards中搜索lung cancer得到的相关基因，组成ground truth基因集，同时计算所选社团的超几何分布检验的p值
% 导入 MalaCards 肺癌基因集合（a,b,c,d,及其总集）为变量 GT_gene，
% 导入 _ENid_BP_GSE19804_Final.mat 中的 gene_id_symbol
% 挑选 Wscore 最大的社团的基因并存为 Cluster_genes。
for i=1:10
    GT_gene_1=cellstr(GT_gene);
    gene_id_symbol_1=cellstr(gene_id_symbol);
    Cluster_genes_1=cellstr(Cluster_genes);
    K=length(GT_gene_1);
    N=length(Cluster_genes_1);
    NK=length(intersect(Cluster_genes_1,GT_gene_1));
    NK_1=length(intersect(gene_id_symbol_1,GT_gene_1));
%     M=length(gene_id_symbol);
    P_value=sum(hygepdf(0:(NK-1),length(gene_id_symbol),K,N)); % Y=hygepdf(X,M,K,N)
end
% for i=1:10
%     GT_gene_1=cellstr(GT_gene);
%     gene_id_symbol_1=cellstr(gene_id_symbol);
%     K=length(intersect(gene_id_symbol_1,GT_gene_1));
%     N=length(Cluster_genes);
%     %
%     P_value=sum(hygepdf(0:(K-1),length(gene_id_symbol),K,N)); % Y=hygepdf(X,M,K,N)
% %     fprintf(fid,'%s\n',['P-value: ' num2str(P_value)]);
% %     
% %     fclose(fid);%关闭文件
% end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/8/24%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 实验三：FLE 评估结果质量
% 原始十个模块的标号以及全网络模块划分结果均保存于
% MyRMOEA__ENid_BP_GSE19804_Final_NewFW_4_C_minWW_ 中，
% 以[0 0 0]为间隔，前面是十个模块（导入为original_Module），后面是总划分(导入为original_CC)。

% 对于RMOEA: Cluster_F=tabulate(original_CC);
            %Selected_Modules=sortrows(Cluster_F,-2);
% 对于SigMod: selected_id=[];
            % for i=1:length(SigMod_genes)
            %     bi_id=cell2mat(cellfun(@(x)SigMod_genes{i,1}==x,gene_id_symbol,'UniformOutput',false));
            %     selected_id=[selected_id;find(bi_id)];
            % end
            % FLE_times=zeros(10,10);
            % for t=1:10
            %     fle = FLE(ones(length(selected_id),1),1,sim_Matrix(selected_id,selected_id),after_PPI_Network(selected_id,selected_id));
            %     FLE_times(t,1:10)=fle;
            % end
% 对于MCODE: selected_id=[];
            % for i=1:length(MCODE_genes)
            %     bi_id=cell2mat(cellfun(@(x)MCODE_genes{i,1}==x,gene_id_symbol,'UniformOutput',false));
            %     selected_id=[selected_id;find(bi_id)];
            % end
            % FLE_times=zeros(10,10);
            % for t=1:10
            %     fle = FLE(original_CC,(1:10),sim_Matrix(selected_id,selected_id),after_PPI_Network(selected_id,selected_id));
            %     FLE_times(t,1:10)=fle;
            % end
% 对于MODA: Cluster_F=tabulate(original_CC);
            %Selected_Modules=sortrows(Cluster_F,-2);
FLE_times=[];
for t=1:10
    % original_CC 需要保证列向量输入。
    % original_Module 需要保证行向量输入。
    fle = FLE(original_CC',original_Module,sim_Matrix,after_PPI_Network);
    FLE_times(t,1:length(fle))=fle;
end
% 复制 FLE_times 于 excel 中


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/8/29%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 实验四：各种基因筛选结束后，使用剩下的基因，用PCC建网后使用PPI筛选边。
% 导入 _ENid_BP_GSE19804_Final.mat
gene_num=length(gene_id_symbol);
Network=zeros(gene_num,gene_num);
% filter_original_Matrix=original_Matrix(logical(node_id),:);
for i=1:gene_num-1
    for j=(1+i):gene_num
        I = corr(original_Matrix(i,:)',original_Matrix(j,:)','type','pearson');
        if after_PPI_Network(i,j)~=0
            Network(i,j)=I;
            Network(j,i)=I;
        end
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/8/23%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 实验五：重新组织基因和模块标号数据，输出，for clusterProfiler GO富集分析
selected_id=[];
for i=1:length(Cluster)
    bi_id=cell2mat(cellfun(@(x)Cluster{i,1}==x,gene_id_symbol,'UniformOutput',false));
    selected_id=[selected_id;find(bi_id)];
end
% 复制 Genes，替换 ENid_BP_GSE18842_Final_gene.txt 的第一列
Genes=sim_genes(selected_id,1);



%% MCODE对比实验
% 导入得分最高的模块的基因为 MCODE_genes[6,1]
for i=1:10
    Accuracy=[];
    Precision=[];
    Recall=[];
    F_measure=[];
    Average_AUC=[];
    % 对于每个社团计算 10 次
    for t=1:10
        Cluster_genes=Cluster(Cluster_id==i);
        selected_id=[];
        for j=1:length(Cluster_genes)
            bi_id=cell2mat(cellfun(@(x)Cluster_genes{j,1}==x,gene_id_symbol,'UniformOutput',false));
            selected_id=[selected_id;find(bi_id)];
        end
        [Indicators,Average_AUC(t,1)] = FiveFold_SVM(original_Matrix(selected_id,:)',label);
        Accuracy(t,1)=Indicators.accuracy;
        Precision(t,1)=Indicators.precision;
        Recall(t,1)=Indicators.recall;
        F_measure(t,1)=Indicators.F_measure;
        name='D:/_Documents/_Papers/Paper_3/Comparison_Methods/MCODE/GSE19188_Result/GSE19188_MCODE_Result_C';
        write_Net=[name num2str(i) '_' 'Indicators.txt'];
        fid=fopen(write_Net,'w');%写入文件路径
        fprintf(fid,'%s\n',['ACC: ' num2str(mean(Accuracy))]);
        fprintf(fid,'%s\n',['Precision: ' num2str(mean(Precision))]);
        fprintf(fid,'%s\n',['Recall: ' num2str(mean(Recall))]);
        fprintf(fid,'%s\n',['F_measure: ' num2str(mean([Indicators.F_measure']))]);
        fprintf(fid,'%s\n',['G-mean: ' num2str(mean([Indicators.G_mean]))]);
        fprintf(fid,'%s\n',['Average_AUC: ' num2str(mean(Average_AUC))]);
    end
    
    %% 实验二：通过MalaCards中搜索lung cancer得到的相关基因，组成ground truth基因集，同时计算所选社团的超几何分布检验的p值
    % 导入 MalaCards 肺癌总基因集合为变量 GT_gene，
    % 导入 _ENid_BP_GSE19804_Final.mat 中的 gene_id_symbol
    % 挑选 Wscore 最大的社团的基因并存为 max_C_gene。
    gene_id_symbol_1=cellstr(gene_id_symbol);
    GT_gene_1=cellstr(GT_gene);
    K=length(intersect(gene_id_symbol_1,GT_gene_1));
    N=length(Cluster_genes);
    P_value=hygepdf(1,length(gene_id_symbol),K,N); % Y=hygepdf(X,M,K,N)
    fprintf(fid,'%s\n',['P-value: ' num2str(P_value)]);
    
    fclose(fid);%关闭文件
end




%% SigMod对比实验
[row,col,w] = find(after_PPI_Network);
final=[gene_id_symbol(col,1) gene_id_symbol(row,1) mat2cell(w,ones(length(w),1),1)];
% 输出数据：前两列为节点的 gene symbol，第三列为边权。
fpath='after_PPI_Network_edge_list.txt';
fid=fopen(fpath,'wt');%写入文件路径
out=final;                            %input_matrix为待输出矩阵
[m,n]=size(out);
for i=1:m
    for j=1:n
        x=out{i,j};
        if j==n
            x=num2str(x);
            fprintf(fid,'%s\n',x);
        else
            fprintf(fid,'%s\t',x);
        end
    end
end
%撤销文件句柄
fclose(fid);

% 生成limma_genes_adjp.txt
% 导入基因的差异表达分析结果为 limma_genes，
% 导入 _ENid_BP_GSE18842_SoftFinal.mat 中的 sim_genes 和 gene_id_symbol
% 在 gene_id 中查找每一个 limma_genes 中的基因，并提取查到的基因的p值
gene_p=zeros(size(gene_id));
for i=1:length(limma_genes)
    for j=1:length(gene_id)
        if cell2mat(limma_genes(i,1))==cell2mat(gene_id(j,1))
            % 将查到的gene的p值保存在gene_p的该基因在gene_id中的位置。
            gene_p(j,1)=cell2mat(limma_genes(i,2));
        end
    end
end
gene_p=[gene_id_symbol num2cell(gene_p)];

fpath='limma_genes_adjp.txt';
fid=fopen(fpath,'wt');%写入文件路径
out=gene_p;                            %input_matrix为待输出矩阵
[m,n]=size(out);
for i=1:m
    for j=1:n
        x=cell2mat(out{i,j});
        if j~=n
            fprintf(fid,'%s\t',x);
        else
            x=num2str(x);
            fprintf(fid,'%s\n',x);
        end
    end
end
%撤销文件句柄
fclose(fid);

% 转入 Rstudio 执行 SigMod，最多可得两个模块（选中的模块及其后续模块）
% 两个模块中的基因分别存储于 selected_genes.tab 和 selected_genes_next.tab
% 导入 selected_genes.tab 为 SigMod_genes，继而执行实验一和实验二
Accuracy=[];
Precision=[];
Recall=[];
F_measure=[];
Average_AUC=[];
% 计算 10 次
for t=1:10
    selected_id=[];
    for i=1:length(SigMod_genes)
        bi_id=cell2mat(cellfun(@(x)SigMod_genes{i,1}==x,gene_id_symbol,'UniformOutput',false));
        selected_id=[selected_id;find(bi_id)];
    end
    [Indicators,Average_AUC(t,1)] = FiveFold_SVM(original_Matrix(selected_id,:)',label);
    Accuracy(t,1)=Indicators.accuracy;
    Precision(t,1)=Indicators.precision;
    Recall(t,1)=Indicators.recall;
    F_measure(t,1)=Indicators.F_measure;
    write_Net='D:/_Documents/_Papers/Paper_3/Comparison_Methods/SigMod_v2/GSE19188/GSE19188_Indicators_SigMod_Result.txt';
    fid=fopen(write_Net,'w');%写入文件路径
    fprintf(fid,'%s\n',['ACC: ' num2str(mean(Accuracy))]);
    fprintf(fid,'%s\n',['Precision: ' num2str(mean(Precision))]);
    fprintf(fid,'%s\n',['Recall: ' num2str(mean(Recall))]);
    fprintf(fid,'%s\n',['F_measure: ' num2str(mean([Indicators.F_measure]))]);
    fprintf(fid,'%s\n',['G-mean: ' num2str(mean([Indicators.G_mean]))]);
    fprintf(fid,'%s\n',['Average_AUC: ' num2str(mean(Average_AUC))]);
    
    % 导入 MalaCards 肺癌总基因集合为变量 GT_gene，
    gene_id_symbol_1=cellstr(gene_id_symbol);
    GT_gene_1=cellstr(GT_gene);
    K=length(intersect(gene_id_symbol_1,GT_gene_1));
    N=length(SigMod_genes);
    P_value=hygepdf(1,length(gene_id_symbol),K,N); % Y=hygepdf(X,M,K,N)
    fprintf(fid,'%s\n',['P-value: ' num2str(P_value)]);
    
    fclose(fid);%关闭文件
end




%% MODA 对比实验
original_Matrix=original_Matrix';

% 共有四个结果：Hierarchical_Result_GSE18842.mat，Hierarchical_Result_GSE19804.mat，Louvain_Result_GSE18842.mat，Louvain_Result_GSE19804.mat
% 分别导入并运行下列代码
% 选择基因个数最大的10个模块
Cluster_F=tabulate(Community);
Selected_Modules=sortrows(Cluster_F,-2);
Selected_Modules(Selected_Modules(:,1)==0,:)=[];
Selected_Modules=Selected_Modules(1:10,1);
for i=1:10
    Accuracy=[];
    Precision=[];
    Recall=[];
    F_measure=[];
    Average_AUC=[];
    for t=1:10
        selected_id=find(Community==Selected_Modules(i,1));
        [Indicators,Average_AUC(t,1)] = FiveFold_SVM(original_Matrix(selected_id,:)',label);
        Accuracy(t,1)=Indicators.accuracy;
        Precision(t,1)=Indicators.precision;
        Recall(t,1)=Indicators.recall;
        F_measure(t,1)=Indicators.F_measure;
        write_Net=['Louvain_Result_GSE19188_Indicators_C' num2str(Selected_Modules(i,1)) '.txt'];
        fid=fopen(write_Net,'w');%写入文件路径
        fprintf(fid,'%s\n',['ACC: ' num2str(mean(Accuracy))]);
        fprintf(fid,'%s\n',['Precision: ' num2str(mean(Precision))]);
        fprintf(fid,'%s\n',['Recall: ' num2str(mean(Recall))]);
        fprintf(fid,'%s\n',['F_measure: ' num2str(mean([Indicators.F_measure]))]);
        fprintf(fid,'%s\n',['G-mean: ' num2str(mean([Indicators.G_mean]))]);
        fprintf(fid,'%s\n',['Average_AUC: ' num2str(mean(Average_AUC))]);
    end
    %% 实验二：通过MalaCards中搜索lung cancer得到的相关基因，组成ground truth基因集，同时计算所选社团的超几何分布检验的p值
    % 导入 MalaCards 肺癌总基因集合为变量 GT_gene，
    % 导入 _ENid_BP_GSE19804_Final.mat 中的 gene_id_symbol
    % 挑选 Wscore 最大的社团的基因并存为 max_C_gene。
    gene_id_symbol_1=cellstr(gene_id_symbol);
    GT_gene_1=cellstr(GT_gene);
    K=length(intersect(gene_id_symbol_1,GT_gene_1));
    N=length(selected_id);
    P_value=hygepdf(1,length(gene_id_symbol),K,N); % Y=hygepdf(X,M,K,N)
    fprintf(fid,'%s\n',['P-value: ' num2str(P_value)]);
    
    fclose(fid);%关闭文件
end