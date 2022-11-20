function W = Wscore(cluster,BinaryAdj)

cnum=max(cluster);
W=zeros(cnum,1);
for i=1:cnum
    Snum=sum(cluster==i);
    Cid=find(cluster==i);
    OS=sum(sum(triu(BinaryAdj(Cid,Cid))));
    BS=0;
    for j=1:length(Cid)
        BS=BS+sum(cluster(BinaryAdj(Cid(1),:))~=i);
    end
%     W(i)=(length(cluster)-Snum)*OS/Snum-BS;
    W(i)=(OS/Snum^2)-(BS/Snum*(length(cluster)-Snum));
% W=sum(W);
end

