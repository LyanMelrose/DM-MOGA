function AdjMatrix=vary_edge_Adj(DATAfile)
fid=fopen(DATAfile,'r');
Df=textscan(fid,'%d %d');
fclose(fid);
Df{1}=Df{1}+1;
Df{2}=Df{2}+1;
Pro=union(Df{1},Df{2});
node_num=length(Pro);
edge_num=length(Df{1});
[bool,edge(:,1)]=ismember(Df{1},Pro);
[bool,edge(:,2)]=ismember(Df{2},Pro);
 AdjMatrix=zeros(node_num);
 for i=1:length(Df{1})
   AdjMatrix(edge(i,2),edge(i,1))=1;
   AdjMatrix(edge(i,1),edge(i,2))=1;
 end
 end
 
