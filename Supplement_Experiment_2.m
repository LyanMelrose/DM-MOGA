%% ���� Main_Process ��ʵ��2�ĳ�����pֵ��������˴���������²���ʵ��
% �ֱ�ȡģ���������ʮ��ȡģ���еĻ�������������֤
% ���� My_RMOEA-��Ȩ ��õ�ʵ����Ϊ ClusterԪ�����飨���򣩣�ͬʱ����ѡȡ��ģ������Ϊ Cluster_id
for i=1:10
% ����ÿ������
        Cluster_genes=Cluster(Cluster_id==i);
        selected_id=[];
        for j=1:length(Cluster_genes)
            bi_id=cell2mat(cellfun(@(x)Cluster_genes{j,1}==x,gene_id_symbol,'UniformOutput',false));
            selected_id=[selected_id;find(bi_id)];
        end
        name='results1/_ENid_BP_GSE19804_Final_C/metrics/_ENid_BP_GSE19804_Final_C';
%         name='D:/_Documents/_Papers/Paper_3/Comparison_Methods/My_RMOEA/results1/_ENid_BP_GSE19188_Final_C/metrics/_ENid_BP_GSE19188_Final_C';
        write_Net=[name num2str(i) '_' 'Pvalue.txt'];
        fid=fopen(write_Net,'w');%д���ļ�·��
    
    %% ʵ�����ͨ��MalaCards������lung cancer�õ�����ػ������ground truth���򼯣�ͬʱ������ѡ���ŵĳ����ηֲ������pֵ
    % ���� MalaCards �ΰ��ܻ��򼯺�Ϊ���� GT_gene��
    % ���� _ENid_BP_GSE19804_Final.mat �е� gene_id_symbol
    % ��ѡ Wscore �������ŵĻ��򲢴�Ϊ max_C_gene��
    gene_id_symbol_1=cellstr(gene_id_symbol);
    GT_gene_1=cellstr(GT_gene);
    Ng=length(GT_gene_1);
    Ngm=length(intersect(cellstr(Cluster_genes),GT_gene_1));
    N=length(Cluster_genes);
    P_value=sum(hygepdf(Ngm:Ng,length(gene_id_symbol),Ng,N)); % Y=hygepdf(X,M,K,N)
    fprintf(fid,'%s\n',['P-value: ' num2str(P_value)]);
    fclose(fid);%�ر��ļ�
end