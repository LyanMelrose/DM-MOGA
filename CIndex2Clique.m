function Clique = CIndex2Clique(Community)
cnum=max(Community);
Clique=cell(cnum,1);
for i=1:cnum
    Clique{i,1}=find(Community==i);
end
end

