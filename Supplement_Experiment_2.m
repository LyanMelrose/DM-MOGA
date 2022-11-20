%% 由于 Main_Process 中实验2的超几何p值计算出现了错误，因此重新补充实验
% 分别取模块个数最大的十个取模块中的基因用来五折验证
% 导入 My_RMOEA-加权 获得的实验结果为 Cluster元胞数组（基因），同时导入选取的模块索引为 Cluster_id
for i=1:10
% 对于每个社团
        Cluster_genes=Cluster(Cluster_id==i);
        selected_id=[];
        for j=1:length(Cluster_genes)
            bi_id=cell2mat(cellfun(@(x)Cluster_genes{j,1}==x,gene_id_symbol,'UniformOutput',false));
            selected_id=[selected_id;find(bi_id)];
        end
        name='results1/_ENid_BP_GSE19804_Final_C/metrics/_ENid_BP_GSE19804_Final_C';
%         name='D:/_Documents/_Papers/Paper_3/Comparison_Methods/My_RMOEA/results1/_ENid_BP_GSE19188_Final_C/metrics/_ENid_BP_GSE19188_Final_C';
        write_Net=[name num2str(i) '_' 'Pvalue.txt'];
        fid=fopen(write_Net,'w');%写入文件路径
    
    %% 实验二：通过MalaCards中搜索lung cancer得到的相关基因，组成ground truth基因集，同时计算所选社团的超几何分布检验的p值
    % 导入 MalaCards 肺癌总基因集合为变量 GT_gene，
    % 导入 _ENid_BP_GSE19804_Final.mat 中的 gene_id_symbol
    % 挑选 Wscore 最大的社团的基因并存为 max_C_gene。
    gene_id_symbol_1=cellstr(gene_id_symbol);
    GT_gene_1=cellstr(GT_gene);
    Ng=length(GT_gene_1);
    Ngm=length(intersect(cellstr(Cluster_genes),GT_gene_1));
    N=length(Cluster_genes);
    P_value=sum(hygepdf(Ngm:Ng,length(gene_id_symbol),Ng,N)); % Y=hygepdf(X,M,K,N)
    fprintf(fid,'%s\n',['P-value: ' num2str(P_value)]);
    fclose(fid);%关闭文件
end