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
All_Run_maxWW=[];   %�洢��ʮ��ʵ�飬ÿ��WWֵ���Ľ����Qֵ��NMIֵ
All_Run_maxWW_Cres=[];   %�洢��ʮ��ʵ�飬ÿ��WWֵ���Ľ����Qֵ��NMIֵ

%% �����ڽӾ����߱�
for mg=1:runtimes
    load(path);
    load('MalaCards_all_connected_genes.mat');
    AdjMatrix=after_PPI_Network;
    BinaryAdj=AdjMatrix~=0;
    degree=uint16(sum(BinaryAdj,1));
%     % ���ڽӾ�����ȥ�����й����ĵ�
%     index=find(degree==0);
%     AdjMatrix(index,:)=[];
%     AdjMatrix(:,index)=[];
    numVar=length(AdjMatrix);
    
    %% ��ȡ�ϲ�����ż���Ӧ����Ϣ
    t3=clock;
    % Ԥ��������
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/6/16%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     [Clique,Matrix,degree1,edgeslist]=local_expansion(AdjMatrix);
    Clique=local_expansion(AdjMatrix);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2020/12/27%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MCODE vertex weighting
    disp('==========VertexWeighting==========')
    Weights=VertexWeighting(BinaryAdj);
    % �õ���ű�ʾ�����Ż��ֽ��
    Pre_Community = Clique2CIndex(Clique);

    Branch_Set=ones(numVar,1);
    % ��ʼ���޼�֦辺�����Ż���
    New_Community=Pre_Community;
    for i=1:length(Clique)
        i
        % ���ڵ�i������
        C_node=cell2mat(Clique(i));
        % ��������д��ڵ���3���ڵ㣬��New_Community��Branch_Set�н���Щ�ڵ�Ķ�Ӧλ�ñ��Ϊ0
        if length(C_node)>=3
%             Branch_Set(C_node,1)=0;
%             New_Community(C_node,1)=0;
%         else
            % ��ȡ��i�����ŵĽڵ㼰��Ȩ�أ��Ը����Ž����޼�,�������Ľڵ���New_Community��Branch_Set�б��Ϊ0
            [New_Clique,Branch_node]=Cut_branch(AdjMatrix,C_node,Weights(C_node),i,length(Clique));
            Branch_Set(Branch_node,1)=0;
            New_Community(Branch_node)=0;
            % ��ȡ�ı����ŵĽڵ㣬�Լ����޼��Ľڵ�
            ChangeID=New_Clique(New_Clique(:,2)~=New_Clique(:,3),1);
            New_Community(ChangeID,1)=New_Clique(New_Clique(:,2)~=New_Clique(:,3),3);
        end
    end
    Cut_Clique=CIndex2Clique(New_Community);
    Cut_Clique(cellfun(@isempty,Cut_Clique))=[];
    % �ҵ����б��޼��Ľڵ㣬�������ڵ�ֱ���Ϊһ���ֲ����ż���Clique�б���
    Branch_node=find(Branch_Set==0);
    Clique_count=length(Cut_Clique);
    for i=1:length(Branch_node)
        Cut_Clique{Clique_count+1}=Branch_node(i);
        Clique_count=Clique_count+1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/6/16%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ����ѹ��ʱ��ô����ߵ�Ȩ�غ���Ҫ��
    Matrix=flod(BinaryAdj,AdjMatrix,Cut_Clique);
    Matrix=single(Matrix);
    
    % edgeslist�а���������ֵ�У����ܵĽ�������ж���
    % һ������ѹ����ʱ������޸�bug��
    % �����������������޸ĺ�����������Ӧ��
    % ��˹ؼ����ж����Ƿ�Ϊ��������ȷ�Ϲ����ľֲ������еĽڵ��Ƿ����û�������ߡ�
    edgeslist = edges_list(Matrix,length(Cut_Clique),Cut_Clique);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/1/4%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    edgeslist([edgeslist(:).n]==0)=[];

    % ָ��'UniformOutput'��ֵΪfalse��cellfun()����Ԫ���������͵Ľ����
    degree1=cell2mat(cellfun(@(S)sum(degree(1,S)),Cut_Clique,'UniformOutput',false));


   %% �ļ����ݱ��漰�����õ�parate��
    root = sprintf('results1/%s/ParetoFront',strNetwork);  %%�����ļ�����ʵ������
    if ~isdir(root) %�ж�·���Ƿ����
        mkdir(root);
    end
    root = sprintf('results1/%s/metrics',strNetwork);
    if ~isdir(root) %�ж�·���Ƿ����
        mkdir(root);
    end
    
    %�����Thirty_Run��ʵ�洢����ÿ��ʵ��Ľ������������ʮ��ʵ��Ľ�����У���ʮ����ʵ����ܽ����All_Run��
    Thirty_Run_maxWW=0;   %�洢����ʵ��Qֵ���Ľ����Qֵ��NMIֵ��RIֵ
    
    fprintf('��%s������\n',num2str(mg));
    
    [C_result,]=TMOEAD(Matrix,Cut_Clique,BinaryAdj,mg,degree1,starttime,edgeslist);
    
    %������ʵ���Thirty_Run������All_Run�С�
    All_Run_maxWW=[All_Run_maxWW;Thirty_Run_maxWW];
    All_Run_maxWW_Cres=[All_Run_maxWW_Cres;C_result];
    
    fprintf('\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/8/16%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,index1]=max(All_Run_maxWW(:,1));    %�ҵ����� Weighted Wscore ����

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/8/16%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ѡ������������10��ģ��
Cluster_F=tabulate(All_Run_maxWW_Cres(index1,:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/8/17%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Selected_Modules=sortrows(Cluster_F,-2);
Selected_Modules(Selected_Modules(:,1)==0,:)=[];
Selected_Modules=Selected_Modules(1:10,1);
% Selected_Modules=Cluster_F(Cluster_F(:,2)>100,1);
%Selected_Modules(Selected_Modules==0)=[];

% ���εõ�ÿ����ѡ�е�ģ��Ļ����������һ��Ϊ gene symbol���ڶ���Ϊģ���ţ���ԭʼ��Ÿ�Ϊ1234˳�ţ�
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
