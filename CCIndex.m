function NeibC = CCIndex( Adj_M )
% Clustering Coefficient of one node


for i=1:length(Adj_M)
    neib=find(Adj_M(i,:)~=0);
    neib_num=length(neib);
%     Nneib=Adj_Table(neib(i),1:len(neib(i)));
    dM=sum(sum(Adj_M(neib,neib)));
    NeibC(i)=dM/neib_num*(neib_num-1);
end
NeibC=sum(NeibC);
end

