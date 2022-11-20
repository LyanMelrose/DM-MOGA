% 导入 sim 矩阵计算后剩余的基因为 sim_genes
% 导入 ENid_BP_GSE18842_sim_genes 2 symbol.txt 的前两列为 gene2symbol
% 导入 ENid_GSE19804_limma_genes_index.csv 的第一列和第二列为差异表达基因集合 diff_gene_id
% 查找每一个 sim_genes 中的基因，并提取查到的基因的p值

gene_p=zeros(size(sim_genes));
% 将 sim_genes 中的 ENTREZ ID 转化为 gene symbol，
gene_id_symbol=cell(size(sim_genes));     
for i=1:length(sim_genes)
    diff_id=find(cell2mat(sim_genes(i,1))==cell2mat(diff_gene_id(:,1)));
    % 将查到的gene的p值保存在gene_p的该基因在sim_genes中的位置。
    if ~isempty(diff_id)
        gene_p(i,1)=cell2mat(diff_gene_id(diff_id,2));
    end
    symbol_id=find(cell2mat(sim_genes(i,1))==cell2mat(gene2symbol(:,1)));
    if ~isempty(symbol_id)
        gene_id_symbol(i,1)=gene2symbol(symbol_id,2);
    end
end

gene_p=[gene_id_symbol num2cell(gene_p)];
fpath='gene_p.txt';
fid=fopen(fpath,'wt');%写入文件路径
out=gene_p;                            %input_matrix为待输出矩阵
[m,n]=size(out);
for i=1:m
    for j=1:n
        x=cell2mat(out(i,j));
        if j==n
            fprintf(fid,'%s\n',x);
        else
            x=num2str(x);
            fprintf(fid,'%s\t',x);
        end
    end
end
%撤销文件句柄
fclose(fid);

