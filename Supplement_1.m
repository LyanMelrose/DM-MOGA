for i=1:10
    GT_gene_1=cellstr(GT_gene);
    gene_id_symbol_1=cellstr(gene_id_symbol);
    Cluster_genes_1=cellstr(Cluster_genes);
    K=length(GT_gene_1);
    N=length(Cluster_genes_1);
    NK=length(intersect(Cluster_genes_1,GT_gene_1));
%     M=length(gene_id_symbol);
    %
    P_value=sum(hygepdf(0:(NK-1),length(gene_id_symbol),K,N)); % Y=hygepdf(X,M,K,N)
end