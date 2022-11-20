function AdjMatrix=flod(BinaryAdj,matrix,clique)

% �õ����нڵ�Ķ�
degree=sum(BinaryAdj,1);
% �õ��ֲ����ŵĸ���
len=length(clique);
% �õ��ڵ����
numVar=length(BinaryAdj);
% ��ʼ��ѹ������ڽӾ���
AdjMatrix=zeros(len,len);
% ��ʼ���ڵ����ű��������
label=zeros(1,numVar);
for i=1:len
    label(clique{i})=i;
end

% �ڵ���ھӵı�ǩ
for i=1:numVar
    D{i}=label(find(BinaryAdj(i,:)==1));
end

% ����ÿ���ֲ�����
for i=1:len
    % ��ȡ�ֲ�����i�еĽڵ�
    a=clique{i};
    % CΪ�ֲ�����i���ھӾֲ�����
    C=[];
    for j=1:length(a)
        C=[C D{a(j)}];
    end
    C=unique(C);
    neighbor=unique(C);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/8/16%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ��i�ڲ��ı���
    AdjMatrix(i,i)=sum(sum(matrix(clique{i},clique{i})))/2;
    %     AdjMatrix(i,i)=sum(sum(BinaryAdj(clique{i},clique{i})))/2;
    
    % ���ھֲ�����i��ÿ���ھ�����a���������֮������ӱ���
    for j=1:length(C)
        a=C(j);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2020/12/23%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        AdjMatrix(a,i)=sum(sum(matrix(clique{i},clique{a})));
%         AdjMatrix(a,i)=sum(sum(BinaryAdj(clique{i},clique{a})));
        AdjMatrix(i,a)=AdjMatrix(a,i);
    end
    AdjMatrix(i,i)=AdjMatrix(i,i)/2;
end

end