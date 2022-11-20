%% ʹ��R��limma���в�������������ͬʱ���� sim_Matrix��
% ���� sim �����Ծ���Ϊ sim_Matrix��
% ���� sim_genes_expression �ľ���Ϊ original_Matrix����������Ϊ sample_id��
% ������Ϊ sim_genes


%% ʹ�ø�˹��MI���������������໥����
% ���� _ENid_BP_GSE18842_SoftFinal.mat �е� original_Matrix��ת�����ļ���./robince-gcmi-f14aae8/matlab
[gene_num,~]=size(original_Matrix);
Network=zeros(gene_num,gene_num);
for i=1:gene_num-1
    for j=(1+i):gene_num
        %I = kernelmi(Matrix(i,:),Matrix(j,:));
        %��֤���������ݵ���Ϊ��������Ϊ����
        I = gcmi_cc(original_Matrix(i,:)',original_Matrix(j,:)');
        Network(i,j)=I;
        Network(j,i)=I;
    end
end

%min-max��׼��
Net_2=(Network-min(Network(:)))/(max(Network(:))-min(Network(:)));
%Net_2Ϊ[0,1]
% Net_2=single(Net);
% Net_2(logical(eye(size(Net_2))))=1;
%Net_1Ϊ[-1��1]
% Net_1=(Net_2-0.5).*2;
% Net_1(logical(eye(size(Net))))=0;
%I = gcmi_cc(x, y)Ҳ������
Network=Net_2;

%% ͨ�� HPRD �� PPI ɸѡ����ֵ���������ı�
% 1. ���� select_Diff_gene_p.m���� sim_genes �� ENTREZ ID ת��Ϊ gene symbol,
% 2. ���� select_HPRD_Edges_from_Adj.m���õ� PPI_Adj��
% 3. �� PPI_Adj ��������һ����õ� Network��
PPI_Adj=logical(PPI_Adj);
Net_2=zeros(size(Network));
Net_2(PPI_Adj)=Network(PPI_Adj);
Network=Net_2;

% ȥ�� Network �ж�Ϊ0�Ľڵ�
non_zero_id=find(sum(Network)~=0);
Network=Network(non_zero_id,non_zero_id);
% �� gene_id_symbol��gene_p��original_Matrix��sim_genes��sim_Matrix ��ȥ����Ϊ 0 �Ľڵ㡣
% non_zero_id=find(sum(PPI_Adj)~=0);
gene_id_symbol=gene_id_symbol(non_zero_id,1);
gene_p=gene_p(non_zero_id,1);
sim_genes=sim_genes(non_zero_id,1);
original_Matrix=original_Matrix(non_zero_id,:);
sim_Matrix=sim_Matrix(non_zero_id,non_zero_id);
% �̶���Network ����Ϊ after_PPI_Network
% after_PPI_Network��gene_id_symbol��gene_p��original_Matrix��sim_genes��sim_Matrix��non_zero_id(node_id)
% ��������Ϊ _ENid_BP_GSE19804_Final.mat


%% ת��Rstudio��ʹ��WGCNAѡȡ��ֵ
% ��������ֵȡΪ1�����Կ���ʡ����һ����


%% ִ��My_RMOEA - ��Ȩ�㷨
% ���� _ENid_BP_GSE19804_Final.mat��ת����My_RMOEA - ��Ȩ���ļ���
my_run


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/8/2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ʵ��һ��ͨ�����۽�����֤��SVM�����������࣬�Զ���������ѡģ��ķ���Ч����
% �ֱ�ȡģ���������ʮ��ȡģ���еĻ�������������֤
% ���� My_RMOEA-��Ȩ ��õ�ʵ����Ϊ ClusterԪ�����飨���򣩣�ͬʱ����ѡȡ��ģ������Ϊ Cluster_id
for i=1:10
    Accuracy=[];
    Precision=[];
    Recall=[];
    F_measure=[];
    Average_AUC=[];
    % ����ÿ�����ż��� 10 ��
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
        fid=fopen(write_Net,'w');%д���ļ�·��
        fprintf(fid,'%s\n',['ACC: ' num2str(mean(Accuracy))]);
        fprintf(fid,'%s\n',['Precision: ' num2str(mean(Precision))]);
        fprintf(fid,'%s\n',['Recall: ' num2str(mean(Recall))]);
        fprintf(fid,'%s\n',['F_measure: ' num2str(mean([Indicators.F_measure']))]);
        fprintf(fid,'%s\n',['G-mean: ' num2str(mean([Indicators.G_mean]))]);
        fprintf(fid,'%s\n',['Average_AUC: ' num2str(mean(Average_AUC))]);
    end
end
    
%% ʵ�����ͨ��MalaCards������lung cancer�õ�����ػ������ground truth���򼯣�ͬʱ������ѡ���ŵĳ����ηֲ������pֵ
% ���� MalaCards �ΰ����򼯺ϣ�a,b,c,d,�����ܼ���Ϊ���� GT_gene��
% ���� _ENid_BP_GSE19804_Final.mat �е� gene_id_symbol
% ��ѡ Wscore �������ŵĻ��򲢴�Ϊ Cluster_genes��
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
% %     fclose(fid);%�ر��ļ�
% end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/8/24%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ʵ������FLE �����������
% ԭʼʮ��ģ��ı���Լ�ȫ����ģ�黮�ֽ����������
% MyRMOEA__ENid_BP_GSE19804_Final_NewFW_4_C_minWW_ �У�
% ��[0 0 0]Ϊ�����ǰ����ʮ��ģ�飨����Ϊoriginal_Module�����������ܻ���(����Ϊoriginal_CC)��

% ����RMOEA: Cluster_F=tabulate(original_CC);
            %Selected_Modules=sortrows(Cluster_F,-2);
% ����SigMod: selected_id=[];
            % for i=1:length(SigMod_genes)
            %     bi_id=cell2mat(cellfun(@(x)SigMod_genes{i,1}==x,gene_id_symbol,'UniformOutput',false));
            %     selected_id=[selected_id;find(bi_id)];
            % end
            % FLE_times=zeros(10,10);
            % for t=1:10
            %     fle = FLE(ones(length(selected_id),1),1,sim_Matrix(selected_id,selected_id),after_PPI_Network(selected_id,selected_id));
            %     FLE_times(t,1:10)=fle;
            % end
% ����MCODE: selected_id=[];
            % for i=1:length(MCODE_genes)
            %     bi_id=cell2mat(cellfun(@(x)MCODE_genes{i,1}==x,gene_id_symbol,'UniformOutput',false));
            %     selected_id=[selected_id;find(bi_id)];
            % end
            % FLE_times=zeros(10,10);
            % for t=1:10
            %     fle = FLE(original_CC,(1:10),sim_Matrix(selected_id,selected_id),after_PPI_Network(selected_id,selected_id));
            %     FLE_times(t,1:10)=fle;
            % end
% ����MODA: Cluster_F=tabulate(original_CC);
            %Selected_Modules=sortrows(Cluster_F,-2);
FLE_times=[];
for t=1:10
    % original_CC ��Ҫ��֤���������롣
    % original_Module ��Ҫ��֤���������롣
    fle = FLE(original_CC',original_Module,sim_Matrix,after_PPI_Network);
    FLE_times(t,1:length(fle))=fle;
end
% ���� FLE_times �� excel ��


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/8/29%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ʵ���ģ����ֻ���ɸѡ������ʹ��ʣ�µĻ�����PCC������ʹ��PPIɸѡ�ߡ�
% ���� _ENid_BP_GSE19804_Final.mat
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
%% ʵ���壺������֯�����ģ�������ݣ������for clusterProfiler GO��������
selected_id=[];
for i=1:length(Cluster)
    bi_id=cell2mat(cellfun(@(x)Cluster{i,1}==x,gene_id_symbol,'UniformOutput',false));
    selected_id=[selected_id;find(bi_id)];
end
% ���� Genes���滻 ENid_BP_GSE18842_Final_gene.txt �ĵ�һ��
Genes=sim_genes(selected_id,1);



%% MCODE�Ա�ʵ��
% ����÷���ߵ�ģ��Ļ���Ϊ MCODE_genes[6,1]
for i=1:10
    Accuracy=[];
    Precision=[];
    Recall=[];
    F_measure=[];
    Average_AUC=[];
    % ����ÿ�����ż��� 10 ��
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
        fid=fopen(write_Net,'w');%д���ļ�·��
        fprintf(fid,'%s\n',['ACC: ' num2str(mean(Accuracy))]);
        fprintf(fid,'%s\n',['Precision: ' num2str(mean(Precision))]);
        fprintf(fid,'%s\n',['Recall: ' num2str(mean(Recall))]);
        fprintf(fid,'%s\n',['F_measure: ' num2str(mean([Indicators.F_measure']))]);
        fprintf(fid,'%s\n',['G-mean: ' num2str(mean([Indicators.G_mean]))]);
        fprintf(fid,'%s\n',['Average_AUC: ' num2str(mean(Average_AUC))]);
    end
    
    %% ʵ�����ͨ��MalaCards������lung cancer�õ�����ػ������ground truth���򼯣�ͬʱ������ѡ���ŵĳ����ηֲ������pֵ
    % ���� MalaCards �ΰ��ܻ��򼯺�Ϊ���� GT_gene��
    % ���� _ENid_BP_GSE19804_Final.mat �е� gene_id_symbol
    % ��ѡ Wscore �������ŵĻ��򲢴�Ϊ max_C_gene��
    gene_id_symbol_1=cellstr(gene_id_symbol);
    GT_gene_1=cellstr(GT_gene);
    K=length(intersect(gene_id_symbol_1,GT_gene_1));
    N=length(Cluster_genes);
    P_value=hygepdf(1,length(gene_id_symbol),K,N); % Y=hygepdf(X,M,K,N)
    fprintf(fid,'%s\n',['P-value: ' num2str(P_value)]);
    
    fclose(fid);%�ر��ļ�
end




%% SigMod�Ա�ʵ��
[row,col,w] = find(after_PPI_Network);
final=[gene_id_symbol(col,1) gene_id_symbol(row,1) mat2cell(w,ones(length(w),1),1)];
% ������ݣ�ǰ����Ϊ�ڵ�� gene symbol��������Ϊ��Ȩ��
fpath='after_PPI_Network_edge_list.txt';
fid=fopen(fpath,'wt');%д���ļ�·��
out=final;                            %input_matrixΪ���������
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
%�����ļ����
fclose(fid);

% ����limma_genes_adjp.txt
% �������Ĳ�����������Ϊ limma_genes��
% ���� _ENid_BP_GSE18842_SoftFinal.mat �е� sim_genes �� gene_id_symbol
% �� gene_id �в���ÿһ�� limma_genes �еĻ��򣬲���ȡ�鵽�Ļ����pֵ
gene_p=zeros(size(gene_id));
for i=1:length(limma_genes)
    for j=1:length(gene_id)
        if cell2mat(limma_genes(i,1))==cell2mat(gene_id(j,1))
            % ���鵽��gene��pֵ������gene_p�ĸû�����gene_id�е�λ�á�
            gene_p(j,1)=cell2mat(limma_genes(i,2));
        end
    end
end
gene_p=[gene_id_symbol num2cell(gene_p)];

fpath='limma_genes_adjp.txt';
fid=fopen(fpath,'wt');%д���ļ�·��
out=gene_p;                            %input_matrixΪ���������
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
%�����ļ����
fclose(fid);

% ת�� Rstudio ִ�� SigMod�����ɵ�����ģ�飨ѡ�е�ģ�鼰�����ģ�飩
% ����ģ���еĻ���ֱ�洢�� selected_genes.tab �� selected_genes_next.tab
% ���� selected_genes.tab Ϊ SigMod_genes���̶�ִ��ʵ��һ��ʵ���
Accuracy=[];
Precision=[];
Recall=[];
F_measure=[];
Average_AUC=[];
% ���� 10 ��
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
    fid=fopen(write_Net,'w');%д���ļ�·��
    fprintf(fid,'%s\n',['ACC: ' num2str(mean(Accuracy))]);
    fprintf(fid,'%s\n',['Precision: ' num2str(mean(Precision))]);
    fprintf(fid,'%s\n',['Recall: ' num2str(mean(Recall))]);
    fprintf(fid,'%s\n',['F_measure: ' num2str(mean([Indicators.F_measure]))]);
    fprintf(fid,'%s\n',['G-mean: ' num2str(mean([Indicators.G_mean]))]);
    fprintf(fid,'%s\n',['Average_AUC: ' num2str(mean(Average_AUC))]);
    
    % ���� MalaCards �ΰ��ܻ��򼯺�Ϊ���� GT_gene��
    gene_id_symbol_1=cellstr(gene_id_symbol);
    GT_gene_1=cellstr(GT_gene);
    K=length(intersect(gene_id_symbol_1,GT_gene_1));
    N=length(SigMod_genes);
    P_value=hygepdf(1,length(gene_id_symbol),K,N); % Y=hygepdf(X,M,K,N)
    fprintf(fid,'%s\n',['P-value: ' num2str(P_value)]);
    
    fclose(fid);%�ر��ļ�
end




%% MODA �Ա�ʵ��
original_Matrix=original_Matrix';

% �����ĸ������Hierarchical_Result_GSE18842.mat��Hierarchical_Result_GSE19804.mat��Louvain_Result_GSE18842.mat��Louvain_Result_GSE19804.mat
% �ֱ��벢�������д���
% ѡ������������10��ģ��
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
        fid=fopen(write_Net,'w');%д���ļ�·��
        fprintf(fid,'%s\n',['ACC: ' num2str(mean(Accuracy))]);
        fprintf(fid,'%s\n',['Precision: ' num2str(mean(Precision))]);
        fprintf(fid,'%s\n',['Recall: ' num2str(mean(Recall))]);
        fprintf(fid,'%s\n',['F_measure: ' num2str(mean([Indicators.F_measure]))]);
        fprintf(fid,'%s\n',['G-mean: ' num2str(mean([Indicators.G_mean]))]);
        fprintf(fid,'%s\n',['Average_AUC: ' num2str(mean(Average_AUC))]);
    end
    %% ʵ�����ͨ��MalaCards������lung cancer�õ�����ػ������ground truth���򼯣�ͬʱ������ѡ���ŵĳ����ηֲ������pֵ
    % ���� MalaCards �ΰ��ܻ��򼯺�Ϊ���� GT_gene��
    % ���� _ENid_BP_GSE19804_Final.mat �е� gene_id_symbol
    % ��ѡ Wscore �������ŵĻ��򲢴�Ϊ max_C_gene��
    gene_id_symbol_1=cellstr(gene_id_symbol);
    GT_gene_1=cellstr(GT_gene);
    K=length(intersect(gene_id_symbol_1,GT_gene_1));
    N=length(selected_id);
    P_value=hygepdf(1,length(gene_id_symbol),K,N); % Y=hygepdf(X,M,K,N)
    fprintf(fid,'%s\n',['P-value: ' num2str(P_value)]);
    
    fclose(fid);%�ر��ļ�
end