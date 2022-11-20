function CIndex = Clique2CIndex(Clique)
% 用于将cell格式的社团表示转化为普通社团标号的表示形式

global numVar

CIndex=zeros(numVar,1);
CIndex_print=1;
for i=1:length(Clique)
    C_node=cell2mat(Clique(i));
    CIndex(C_node,1)=i;
end
end

