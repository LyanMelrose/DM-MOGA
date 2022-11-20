function SNPTable=CutCommunities(Commu_Adj,SNPTable,CommunitiesTag,seed)
%%

Neighbors=intersect(find(Commu_Adj(seed,:)~=0),find(SNPTable(:,3)==0));
SNPTable(Neighbors,3)=CommunitiesTag;
for i=1:length(Neighbors)
    SNPTable=CutCommunities(Commu_Adj,SNPTable,CommunitiesTag,Neighbors(i));
end