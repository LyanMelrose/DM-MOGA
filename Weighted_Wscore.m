function W = Weighted_Wscore(cluster,Network)
cnum=max(cluster);
W=zeros(cnum,1);
for i=1:cnum
    Snum=sum(cluster==i);
    Cid=logical(cluster==i);
    OS=sum(sum(triu(Network(Cid,Cid))));
    BS=0;
    BS=sum(sum(Network(Cid,~Cid)));
    W(i,1)=(length(cluster)-Snum)*OS/Snum-BS;
%     W(i,1)=OS;
end
% W=max(W);
W=sum(W);
end

