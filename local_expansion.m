function clique=local_expansion(Matrix)

numVar=size(Matrix,1);
 % over_node���ؾֲ����ż��ص��Ľڵ㣬clique���������ҵ��ľֲ�����
[over_node,clique]=localexpansion(Matrix,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2020/6/4%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
over_node=uint16(over_node);
% find_overlap_node()���������ù�
t1=clock;
[clique]=find_overlap_node(clique,over_node,Matrix);


% Q��ͳ�����б����ֵĽڵ�
T=length(clique);
Q=[];
for i=1:length(clique)
    Q=[Q uint16(clique{i})];
end
% ��δ������ֲ����ŵĵ㲹�ھֲ���������clique�����
A=1:length(Matrix);
A=setdiff(A,Q);
for i=1:length(A)
    clique{T+1}=A(i);
    T=T+1;
end

clique(cellfun(@isempty,clique))=[];
t=length(clique);
% ��Ԫ����������
% ΪʲôҪ����
[~,rank] = sort(cell2mat(cellfun(@(s)s(1),clique,'UniformOutput',false)));
clique= clique(rank);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2020/12/24%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% degree=sum(Matrix,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/1/2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ���������ԣ���֤cliqueȷʵ���ص����ֵĽڵ㡣
% count=zeros(numVar,1);
% for i=1:length(clique)
%     cc=cell2mat(clique(i));
%     count(cc,1)=count(cc,1)+1;
% end

t1=clock;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/6/16%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ����ʡ�Դ�����ԭ������
% Matrix=flod(Matrix,clique);
% Matrix=single(Matrix);
% edgeslist = edges_list(Matrix,t,clique);
% degree1=cell2mat(cellfun(@(S)sum(degree(1,S)),clique,'UniformOutput',false));

end
