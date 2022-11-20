function  [SNPTable,Branch_node] = Cut_branch(Adj,C_node,Weights,C_Index,Cnum)

%%
[A,B]=sort(Weights,'descend');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/8/28%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Avg_Cnode=mean(sum(Adj(C_node,C_node)));
% 删掉权重小于等于2，且
if sum(A>2)~=0
    Index_a=C_node(B(A>2));
    Index_b=C_node(sum(Adj(C_node,C_node))>=Avg_Cnode);
    Index=intersect(Index_a,Index_b);
else
    Index=C_node;
end
% % 删掉权重小于等于2的节点
% if sum(A>2)~=0
%     Index=C_node(B(A>2));
% else
%     Index=C_node;
% end

% 提取社团内部专属邻接矩阵
Commu_Adj=Adj(Index,Index);
% 节点的标号
SNPTable(:,1)=Index;
% 节点原本所属的社团
SNPTable(:,2)=ones(size(Index))*C_Index;
% 节点目前所属的社团
SNPTable(:,3)=0;

CommunitiesTag=Cnum;%注意多出来的局部社团怎么办

% 当还有节点未被划分新社团时，进行如下循环
while ~isempty(find(SNPTable(:,3)==0, 1))
    candidate_seed=SNPTable(:,3)==0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/3/2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 改为内部链接度
    % 在SNPTable中序号为seed的节点作为种子节点，计算候选种子节点的度
%     candidate_seed_degree=sum(Commu_Adj(:,candidate_seed));
    % 如果仅剩一个未划分新社团的节点，因此被成为候选种子节点，且其内部链接度为0，直接跳出循环，该节点被修剪。
%     if all(candidate_seed_degree)==0
%         break;
%     end

    %% 改为整体度（不分社团内部还是外部）
    candidate_seed_degree=sum(Adj(:,SNPTable(candidate_seed,1)));
    odd_degree=zeros(size(Index));
    odd_degree(candidate_seed)=candidate_seed_degree;

    % 选择度最大的候选节点作为种子节点
    [~,seed]=max(odd_degree);
    
    CommunitiesTag=CommunitiesTag+1;
    SNPTable(seed,3)=CommunitiesTag;
    % 从种子节点递归地将其待划分新社团的邻居归入社团CommunitiesTag中
    SNPTable=CutCommunities(Commu_Adj,SNPTable,CommunitiesTag,seed);
end

% 统计新生成的社团及每个新社团的大小
New_Commu=unique(SNPTable(:,3));
count=zeros(size(New_Commu));
for i=1:length(New_Commu)
    count(i,1)=sum(SNPTable(:,3)==New_Commu(i,1));
end
% 选择最大的新社团，将其设为原社团，
[~,ori_I]=max(count);
ori_I=New_Commu(ori_I,1);
SNPTable(SNPTable(:,3)==ori_I,3)=C_Index;
New_Commu(New_Commu==ori_I)=C_Index;
% 如果新社团中存在大小低于3的社团，将其作为被删除的社团，其中的节点在SNPTable的第三列被标记为0。
if sum(count<3)~=0
    Delete_New_Commu=New_Commu(count<3);
    for i=1:length(Delete_New_Commu)
        SNPTable(SNPTable(:,3)==Delete_New_Commu(i,1),3)=0;
    end
end

% setdiff()返回在C_node中，但不在Index中的元素，即在函数最开始因为权重小于等于二被修剪的节点
Branch_node=setdiff(C_node,Index);
% 再补充SNPTable的第三列为0的那些节点，共同组成了被修剪的节点集合。
Branch_node=[Branch_node SNPTable(SNPTable(:,3)==0,1)'];
end