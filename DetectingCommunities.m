function Community_divide=DetectingCommunities(Weights,VertexTable,AdjacentMatrix)
%%
[A,B]=sort(Weights,'descend');
if sum(A>2)~=0
    Index=B(A>2);
else
    Index=B;
end

SNPs=VertexTable(Index,1);
Communities=AdjacentMatrix(Index,Index);
SNPTable(:,1)=SNPs;
SNPTable(:,2)=1:length(SNPs);
SNPTable(:,3)=0;

CommunitiesTag=0;
%找到第一个未被划分社团的节点的标号，作为种子节点
while ~isempty(find(SNPTable(:,3)==0, 1))
    seed=find(SNPTable(:,3)==0, 1);
    CommunitiesTag=CommunitiesTag+1;
    SNPTable(seed,3)=CommunitiesTag;
    SNPTable=CutCommunities(Communities,SNPTable,CommunitiesTag,seed);
end

Community_divide=SNPTable(:,[1,3]);

%按社团序号给节点排序
SNPTable=SNPTable(:,[3,1,2]);
[~,B]=sort(SNPTable(:,1),'ascend');
SNPTable=SNPTable(B,:);