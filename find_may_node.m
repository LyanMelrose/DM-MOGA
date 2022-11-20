function fault_node=find_may_node(Node3,AdjMatrix,V,numVar)

g=ones(1,numVar)                                 %-----记录散点
for i=1:V
    if length(Node3(i).e)>=3
        g(Node3(i).e)=0;
    end
    
end

fault_node=[];                                 %求在外面的散点数量/在团内的数量>0.5
for u=1:V
if length(Node3(u).e)>=3
for j=1:length(Node3(u).e)
a=find(AdjMatrix(Node3(u).e(j),:)==1);
c= intersect(a,Node3(u).e);
b=setdiff(a,Node3(u).e);
b=b(find(g(b)==1));

if length(b)/length(c)>=1
fault_node=[fault_node Node3(u).e(j)];
end
end
end
end


end