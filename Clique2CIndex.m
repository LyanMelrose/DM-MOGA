function CIndex = Clique2CIndex(Clique)
% ���ڽ�cell��ʽ�����ű�ʾת��Ϊ��ͨ���ű�ŵı�ʾ��ʽ

global numVar

CIndex=zeros(numVar,1);
CIndex_print=1;
for i=1:length(Clique)
    C_node=cell2mat(Clique(i));
    CIndex(C_node,1)=i;
end
end

