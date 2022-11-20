% ���� sim ��������ʣ��Ļ���Ϊ sim_genes
% ���� ENid_BP_GSE18842_sim_genes 2 symbol.txt ��ǰ����Ϊ gene2symbol
% ���� ENid_GSE19804_limma_genes_index.csv �ĵ�һ�к͵ڶ���Ϊ��������򼯺� diff_gene_id
% ����ÿһ�� sim_genes �еĻ��򣬲���ȡ�鵽�Ļ����pֵ

gene_p=zeros(size(sim_genes));
% �� sim_genes �е� ENTREZ ID ת��Ϊ gene symbol��
gene_id_symbol=cell(size(sim_genes));     
for i=1:length(sim_genes)
    diff_id=find(cell2mat(sim_genes(i,1))==cell2mat(diff_gene_id(:,1)));
    % ���鵽��gene��pֵ������gene_p�ĸû�����sim_genes�е�λ�á�
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
fid=fopen(fpath,'wt');%д���ļ�·��
out=gene_p;                            %input_matrixΪ���������
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
%�����ļ����
fclose(fid);

