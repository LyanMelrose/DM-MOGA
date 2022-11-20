function C_Result=TMOEAD(Matrix,Clique,BinaryAdj,mg,degree1,starttime,edgeslist)
%% �����벿��
% ��������������������������أ����ڵ�����ģ����������Ŵ��㷨����ҪĿ�ģ����ߴ�������Ƿֿ��ġ�

global popsize numVar  max_gen  M niche V pc pm
global strNetwork weights neighbors
global Thirty_Run_maxWW   Adj_mat   t3
global original_Matrix    label    gene_id_symbol    GT_gene    AdjMatrix

%% ����ǰ�ڴ���--kk��ʾĿ�꺯������--ideal��ʾ�����--weight��neighbours��ʾȨ�غ��ھ�--Chromosomes��ʾ��Ⱥ
kk=2;
popsize=50;
niche=40;


V=size(Matrix,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/8/2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idealp = Inf*ones(1,M);
% idealp = -Inf*ones(1,M);

% ��ʼ����Ⱥ�и����Ȩ�غ��ھӣ�Ӧ���Ƕ�Ŀ��ֽ��һ����
[weights neighbors] = init_weight(popsize, niche);
popsize = size(weights,1);
index=[];

% ��ʼ����Ⱥ
chromosomes = initialize_variables(popsize,edgeslist,M,degree1,Matrix,Clique);
% chromosomes = initialize_variables(popsize,edgeslist,M,ll,degree1,Matrix,Clique);



%% ��������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/7/30%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stop=1;

Save_chromosome_sub_FF=zeros(max_gen,2);
% ��ÿһ������ѡģ��ĸ�������ָ��ֵ���ηֱ�洢Ϊ��������
% �����ѡģ���������10������������
% ���ʼѡ��ģ��ĸ���Ϊ0
CNUM=0;
Module_ACC=[];
Module_AUC=[];
Module_P=[];
All_Upper_Idx=zeros(max_gen,7);
for Gene = 1 : max_gen
    t1=clock;
    % max_gen=100����˽��� Gene=30 �� 60 ʱִ�н��������е�������������
%     if mod(Gene,30)==0&&Gene<=60
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/7/30%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %         stop=stop+1;
%         
%         % [ Node3,chromosomes, degree1, edgeslist,V,fault_node]= Fault_tolerance(Node3,chromosomes, degree1, edgeslist,V,fault_node,AdjMatrix);            %%����ش����������ĵ㣬����adj_mat��edgeslist��degree1��V
%         Draw(chromosomes,V);
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/7/30%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         [ CLique,chromosomes,idealp,degree1,edgeslist,V,Matrix] = change_evolution(chromosomes,V,Clique,BinaryAdj,degree1,edgeslist,Matrix,AdjMatrix);
%         %         [ Clique,chromosomes, degree1, edgeslist,V,Matrix] = change_evolution(chromosomes,V,Clique,BinaryAdj,ll,degree1,edgeslist,Matrix);
%         Clique=CLique;
%     end
    
    % ������һ�����塣
    % ���ڵ�ǰ��Ⱥ�е�ÿ������
    for i=1:popsize
        % ��ȡ�� i ��������ھ�
        i_neighbor_index = neighbors(i,:);
        % ��ȡ���ھӵ�Ȩ��
        i_neighbor_weight = weights(i_neighbor_index,:);
        % ��ȡ���ھӵ�Ⱦɫ��
        i_neighbor_chromosome = chromosomes(i_neighbor_index,:);
        % ��ȡ���ھӵ���Ӧ�Ⱥ���ֵ
        old_objective = i_neighbor_chromosome(:,V+1:V+M);
        
        % �Ŵ��������ӵĲ�������
        % �����޸���Ⱦɫ�壬ȴû���޸�clique,�Ƿ�Ӧ����change_evolution()�ķ���ֵ��ͬ��
        % ����[ Clique,f,idealp,degree1,edgeslist,V,Matrix]
        % ��̫�ԣ���Ȼ������ģ�黮�ֵ�ת�䣬��δ��������������
        %% ͨ�����桢�����Ŵ��㷨���ӻ�������µ�Ⱦɫ��
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/8/2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        new_chromosomes = genetic_op(i_neighbor_chromosome,chromosomes,edgeslist,pc,pm,V,M,degree1,Matrix,Clique);
        %         new_chromosomes = genetic_op(i_neighbor_chromosome,chromosomes,edgeslist,pc,pm,V,M,ll,degree1,Matrix,Clique);
        
        % ����ȫ����С��Ӧ�Ⱥ���ֵ��update idealp
        idealp = min(min(new_chromosomes(:,V+1:V+M)),idealp);


        % update the neighbours
        % ѡ����ڵ�ǰ��������õĽ�
        chromosome_sub = sort_sub(weights(i,:),new_chromosomes,idealp,V,M);
        Save_chromosome_sub_FF(Gene,:)=chromosome_sub(V+1:V+M);
        
        %-------------------
        new_subobj = subobjective_te(i_neighbor_weight, chromosome_sub(V+1:V+M), idealp);
        old_subobj = subobjective_te(i_neighbor_weight, old_objective, idealp);
        jj=0;
        iposition = find(new_subobj<old_subobj);
        jj=length(iposition);
        if jj>kk         %%��ȡ2���������
            kk1=floor(jj*rand(kk,1))+1;%%rand��2��1�����������2��һ�������
            iposition1=iposition(kk1);
        else
            iposition1=iposition;
        end
        chromosomes(i_neighbor_index(iposition1)',:) = chromosome_sub(ones(length(iposition1),1),:);
        Problem='��������';
        M=2;
        clc;
        fprintf('NSGA-II,��%2s��,%5s,%2sά,�����%4s%%,��ʱ%5s��\n',num2str(1),Problem,num2str(M),num2str(roundn(Gene/max_gen*100,-1)),num2str(roundn(toc,-2)));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/8/10%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ������ѭ������ø������Ϊ��ѹ�����Ż��ֽ����
    D_chromosome_sub=decode1(chromosome_sub(1,1:V),Clique);
    % ������������Ӧ�Ⱥ���ֵ��pֵ�ͷ����ȵĹ�ϵ
    C_Frequency=tabulate(D_chromosome_sub);
    % ѡ��ģ���С����ǰʮ��ģ��
    Selected_Modules=sortrows(C_Frequency,-2);
    Selected_Modules(Selected_Modules(:,1)==0,:)=[];
    Selected_Modules=Selected_Modules(1:10,1);
    % ������ѡģ��ĸ���
    CNUM=10;
    % ����ÿ����ѡ�е�ģ�飬�ֱ����۲���¼ָ��ֵ
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
% ��ͼ,(1,1)ΪACC��(1,2)ΪF1��(2,1)ΪAUC��(2,2)ΪPP
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


%% ������������parate���ϵ�
times=etime(clock,t3);
Population=chromosomes;
%����õ� pareto ǰ��
ParetoFront=unique(Population,'rows');
path = sprintf('results1/%s/ParetoFront/MODPSO1_%s_PF%d.txt',strNetwork,strNetwork,mg);  %%��Ⱥ������ǰ������
%% save_metrics
%ParetoFront1��Ϊ���Ż��ֽ��������ÿ�д���һ��Pareto���Ž⣬ÿ�д���һ���ڵ����������š�
WW=zeros(size(ParetoFront,1),1);
for i=1:size(ParetoFront,1)
    ParetoFront1(i,1:numVar) = decode1(ParetoFront(i,1:V),Clique);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/8/31%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    WW(i,1) = Weighted_Wscore(ParetoFront1(i,:),AdjMatrix);%%����ÿ������ļ�Ȩ W ֵ
end

% metrics=ParetoFront;
[~,index1]=max(WW);
fprintf('Weighted-Wֵ���ʱ���������Ŀ��%d\n',length(unique(ParetoFront1(index1,:))));
C_Result=ParetoFront1(index1,:);

disp(strNetwork);
fprintf('��%d��ʵ����\n',mg);

savedata1(path,ParetoFront);
Thirty_Run_maxWW=WW(index1,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/8/25%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path = sprintf('results1/%s/metrics/%s_All_Upper_Idx_%d.txt',strNetwork,strNetwork,mg);
% ���� All_Upper_Idx
savedata1(path,All_Upper_Idx);
% path = sprintf('results1/%s/metrics/%s_PF_FitnessValue%d.txt',strNetwork,strNetwork,mg);
% savedata1(path,metrics(:,V+1:V+M));


%Ӧ����ʵ�ֵ����һ���������ԣ�ִ����ɺ�������������Thirty_Run_maxWW��path�б�������ݻᱻ�滻������������WW����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/8/16%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
last_operate(strNetwork,BinaryAdj,mg,degree1,Clique);
% last_operate(A,strNetwork,BinaryAdj,mg,degree1,Clique);

total_time = etime(clock, starttime);
path_time = sprintf('results1/%s/metrics/total_time_%d.txt',strNetwork,mg);
savedata1(path_time,total_time);

fprintf('total time used %u\n', etime(clock, starttime));
end