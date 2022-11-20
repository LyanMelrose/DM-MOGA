function  [SNPTable,Branch_node] = Cut_branch(Adj,C_node,Weights,C_Index,Cnum)

%%
[A,B]=sort(Weights,'descend');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/8/28%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Avg_Cnode=mean(sum(Adj(C_node,C_node)));
% ɾ��Ȩ��С�ڵ���2����
if sum(A>2)~=0
    Index_a=C_node(B(A>2));
    Index_b=C_node(sum(Adj(C_node,C_node))>=Avg_Cnode);
    Index=intersect(Index_a,Index_b);
else
    Index=C_node;
end
% % ɾ��Ȩ��С�ڵ���2�Ľڵ�
% if sum(A>2)~=0
%     Index=C_node(B(A>2));
% else
%     Index=C_node;
% end

% ��ȡ�����ڲ�ר���ڽӾ���
Commu_Adj=Adj(Index,Index);
% �ڵ�ı��
SNPTable(:,1)=Index;
% �ڵ�ԭ������������
SNPTable(:,2)=ones(size(Index))*C_Index;
% �ڵ�Ŀǰ����������
SNPTable(:,3)=0;

CommunitiesTag=Cnum;%ע�������ľֲ�������ô��

% �����нڵ�δ������������ʱ����������ѭ��
while ~isempty(find(SNPTable(:,3)==0, 1))
    candidate_seed=SNPTable(:,3)==0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/3/2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% ��Ϊ�ڲ����Ӷ�
    % ��SNPTable�����Ϊseed�Ľڵ���Ϊ���ӽڵ㣬�����ѡ���ӽڵ�Ķ�
%     candidate_seed_degree=sum(Commu_Adj(:,candidate_seed));
    % �����ʣһ��δ���������ŵĽڵ㣬��˱���Ϊ��ѡ���ӽڵ㣬�����ڲ����Ӷ�Ϊ0��ֱ������ѭ�����ýڵ㱻�޼���
%     if all(candidate_seed_degree)==0
%         break;
%     end

    %% ��Ϊ����ȣ����������ڲ������ⲿ��
    candidate_seed_degree=sum(Adj(:,SNPTable(candidate_seed,1)));
    odd_degree=zeros(size(Index));
    odd_degree(candidate_seed)=candidate_seed_degree;

    % ѡ������ĺ�ѡ�ڵ���Ϊ���ӽڵ�
    [~,seed]=max(odd_degree);
    
    CommunitiesTag=CommunitiesTag+1;
    SNPTable(seed,3)=CommunitiesTag;
    % �����ӽڵ�ݹ�ؽ�������������ŵ��ھӹ�������CommunitiesTag��
    SNPTable=CutCommunities(Commu_Adj,SNPTable,CommunitiesTag,seed);
end

% ͳ�������ɵ����ż�ÿ�������ŵĴ�С
New_Commu=unique(SNPTable(:,3));
count=zeros(size(New_Commu));
for i=1:length(New_Commu)
    count(i,1)=sum(SNPTable(:,3)==New_Commu(i,1));
end
% ѡ�����������ţ�������Ϊԭ���ţ�
[~,ori_I]=max(count);
ori_I=New_Commu(ori_I,1);
SNPTable(SNPTable(:,3)==ori_I,3)=C_Index;
New_Commu(New_Commu==ori_I)=C_Index;
% ����������д��ڴ�С����3�����ţ�������Ϊ��ɾ�������ţ����еĽڵ���SNPTable�ĵ����б����Ϊ0��
if sum(count<3)~=0
    Delete_New_Commu=New_Commu(count<3);
    for i=1:length(Delete_New_Commu)
        SNPTable(SNPTable(:,3)==Delete_New_Commu(i,1),3)=0;
    end
end

% setdiff()������C_node�У�������Index�е�Ԫ�أ����ں����ʼ��ΪȨ��С�ڵ��ڶ����޼��Ľڵ�
Branch_node=setdiff(C_node,Index);
% �ٲ���SNPTable�ĵ�����Ϊ0����Щ�ڵ㣬��ͬ����˱��޼��Ľڵ㼯�ϡ�
Branch_node=[Branch_node SNPTable(SNPTable(:,3)==0,1)'];
end