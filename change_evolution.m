function [ Clique,f,idealp,degree1,edgeslist,V,Matrix]=change_evolution(chromosomes,V,CLique,BinaryAdj,degree1,edgeslist,matrix,ori_matrix)
% function [ Clique,f, degree1, edgeslist,V,Matrix]=change_evolution(chromosomes,V,CLique,BinaryAdj,ll,degree1,edgeslist,matrix)
global numVar
% g ΪĿǰ�����нڵ�ĸ���������������ľֲ����ţ�
g=V;      %������������޷��ϲ������
% Ⱥ����Ⱦɫ��ĸ���
popsize=size(chromosomes,1);
% Ŀǰ�����нڵ�ĸ���������������ľֲ����ţ�

% ȡȺ��������Ⱦɫ���������Ӧ�Ⱥ���ֵ
BB=chromosomes(:,[V+1 V+2]);
% ������Ⱦɫ����з�֧������
B=P_sort(BB,'all');
% �õ� pareto ��һǰ�صĳ�Ա����������������Ϊ��Ӣ����
B=find(B==1);
ParetoFront=chromosomes(:,1:V);% ��ȡ����Ⱦɫ������ݣ���������Ӧ�Ⱥ���ֵ��
t1=clock;
%%%�Ը������ڵĵ�����𵴣��ҳ����ڱ仯�ĵ��ܶ�Q�����������á�
% ����Paretoǰ���е����и���
for i=1:size(ParetoFront,1)
    % �������Ϊ����ģ����
    ParetoFront1(i,:) = decode1(ParetoFront(i,1:V),CLique);
    C_num(i,1)=max(ParetoFront1(i,:));   %%ͳ�����и����ģ�黮�ֽ����ģ��ĸ���
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/8/31%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    WW(i,1) = Weighted_Wscore(ParetoFront1(i,:),ori_matrix);%%����ÿ������ļ�Ȩ W ֵ
end


% �ҵ���Ӣ�����л��ִ���Ľڵ�
% �����޸ĺ�ľ�Ӣ����ģ�黮�ֽ�� ParetoFront1(B,:),�Լ�����ڵ������ remove��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/8/16%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ParetoFront1(B,:),remove]=find_error(ParetoFront1(B,:),ori_matrix,CLique,WW);
% [ParetoFront1(B,:),remove]=find_error(ParetoFront1(B,:),matrix,CLique,WW);
%%%��ǩ�ֲ���ɢ----------------
Clique=find_merge(WW,C_num,chromosomes,CLique);
% �����µľֲ����š�
Clique=find_clique(Clique,remove,numVar);
Clique=Clique';
% �õ����оֲ����Ÿ��º󣬾ֲ����ŵĸ��� V
V=length(Clique);
% ����ֲ����Ÿ��������ı�
if V~=g
    % �����ڽӾ�����������
    Matrix=find_Matrix(ori_matrix,Clique);
    % ���±߱���Ϣ�ͽڵ����Ϣ
    edgeslist = edges_list(Matrix,V,Clique);
    degree=sum(BinaryAdj,1);
    degree1=cell2mat(cellfun(@(S)sum(degree(S)),Clique,'UniformOutput',false));
    
    %%%�̳�ǰ��30����ģ�黮����Ϣ������Ȼ�����������������������µľֲ����ţ������нڵ�����ģ�鱣�ֲ���%%%%
    % ����ÿ������
    for i=1:popsize
        % ��ǰ30���и��� i ����ÿ���ֲ����ŵ�ģ�����������þֲ����ŵ����нڵ����ģ������
        for j=1:length(Clique)
            ParetoFront1(i,Clique{j})=ParetoFront1(i,Clique{j}(1));
        end
        % inheritΪ�̳к�������ǰ���޸ĺ�� ParetoFront1 ���룬�����ص������Ǳ����Ľ������chromosomes
        f(i,1:V)=inherit(ParetoFront1(i,:),Clique);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/7/30%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        f(i,V + 1: V+2) = evaluate_objective(f(i,1:V),degree1,Matrix,Clique);
        %         f(i,V + 1: V+2) = evaluate_objective(f(i,1:V),ll,degree1,Matrix,Clique);
    end
    % ������Ӧ�Ⱥ���ȫ����Сֵ
    idealp = min(f(:,V+1:V+2));
else      % ����ֲ����Ÿ���δ�����ı�
    % ������chromosomes��CLique���Լ�matrix
    f=chromosomes;
    idealp = min(f(:,V+1:V+2));
    Clique=CLique;
    Matrix=matrix;
end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ParetoFront1,remove]=find_error(ParetoFront1,AdjMatrix,CLique,W)
M=[];
% V ΪĿǰ�����нڵ�ĸ���������������ľֲ����ţ�
V=length(CLique);
% numVar Ϊδ���������еĽڵ����
numVar=size(AdjMatrix,1);

% ����ÿ������
for in=1:size(ParetoFront1,1)
    % ��ʼ�� change_node
    change_node{in}=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/8/16%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % label����ȡ�� in �������ģ�黮�ֽ��
    label= ParetoFront1(in,1:numVar);
%     label= ParetoFront1(in,1:numVar);
    
    %% %% ����2�ַ�ʽ-����Weighted Wscore label = ParetoFront1(in,1:numVar);
    % ����ÿ���ڵ㣨�����ֲ����ţ�
    for i=1:V
        % ����ֲ������а������������ڵ�
        if(length(CLique{i})>=3)
            % ���ھֲ������е�ÿ���ڵ�
            for j=1:length(CLique{i})
                % �ֲ������б����������ڵ����������ȡ�� j ���ڵ������ m
                m=CLique{i}(j);
                % ��ȡ�� j ���ڵ�����ģ��ı�� k
                k=label(m);
                % index_last���ҵ�ģ�� k �����нڵ������
                index_last=find(label==k);
                % ���ڵ� m �Ƿ�����ھɵľֲ�������
                A=setdiff(index_last,m);
                % neighbors����ȡ�ڵ� m ���ھ�
                neighbors=find(AdjMatrix(m,:)~=0);
                % Max()Ϊ�Զ��庯��������Ѱ�Ұ����ھӸ�������ģ���id
                com_max=Max(neighbors,label);
                % com_max ���ؽڵ� m ���ھ�����ģ�飬��ģ��������ھӵĸ�����ࡣ
                % �������ʹ�ýڵ� m Ҳ������һģ�顣
                ParetoFront1(in,m)=com_max;
                % ����ڵ� m ������ģ�鷢���˸ı�
                if com_max~=k
                    % �� m ���Ϊ�� in ���ֲ����ŵ� change node�������档
                    change_node{in}=[change_node{in} m];
                end
            end
        end
    end
    % ParetoFront1 �Ѿ��ڵ� 101 �н����˸��£������µ� Weighted Wscore
    WW(in) = Weighted_Wscore(ParetoFront1(in,1:numVar),AdjMatrix);
    % ���change node�еĽڵ��޸ĺ��ܹ���ø���� Weighted Wscore ֵ����ȡ�������
    if WW(in)<W(in)                                                       %������
        M=[M in];
    end
end
error_node=[];
% ��M�еĸ����change node ���أ���Ϊ���ִ���Ľڵ㣬�Լ���������ֵ��
for i=1:length(M)
    error_node=[error_node change_node{M(i)}];
end
remove=unique(error_node);

end


function Clique=find_merge(W,C_num,chromosomes,CLique)
%��Ҫ����ָ��Խ��Խ��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/8/16%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f(:,1)=W;
% f(:,1)=-Q;

f(:,2)=-C_num;
V=length(CLique);
FrontValue = P_sort(f,'all');   %front��
B=find(FrontValue==1);
K=chromosomes(B,1:V);
[a,~]=size(K);
MM=[];
for i=1:a
    A=decode(K(i,:));
    MM=[MM;A];
end
visited=zeros(1,V);
times=1;
for i=1:V
    if visited(i)==0
        distance_i=zeros(1,V);
        for j=i+1:V
            if visited(j)==0
                % ͳ�ƾֲ�����i��j��ͬ��λ�������aӦ���Ǹ��������
                A=MM(:,i)-MM(:,j);
                index_0=find(A==0);
                distance_i(j)=length(index_0)/a;
            end
        end
        % �ҵ��ڴ󲿷ָ�������ͬ�ľֲ�����
        F=find(distance_i>0.9);
        visited(F)=1;
        % times��ʲô������
        Clique{times}=cell2mat(CLique([F i]));
        if length(Clique{times})>0
            times=times+1;
        end
    end
end
end



function Clique=find_clique(Clique,erase_node,numVar)
t=length(Clique);
for i=1:t
    Clique{i}=setdiff(Clique{i},erase_node);
end
R=[];
for i=1:t
    R=[R Clique{i}'];
end
R=setdiff(1:numVar,R);
for i=1:length(R)
    t=t+1;
    Clique{t}=R(i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/8/15%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Clique=Clique(find(cell2mat(cellfun(@(S)length(S),Clique,'UniformOutput',false))~=0));  %%�޳�����Ϊ0�����ţ�

[~,rank] = sort(cell2mat(cellfun(@(s)s(1),Clique,'UniformOutput',false)));                %��Ԫ����������
Clique= Clique(rank);
end

function Matrix=find_Matrix(Matrix,Clique)
numVar=size(Matrix,1);
% degree=sum(Matrix,1);
D=1:numVar;
for i=1:length(Clique)
    C=Clique{i};
    vertex_min=find(D==C(1));
    while length(C)>1
        j=2;
        vertex_max=find(D==C(2));
        if length(vertex_max)==0 %�����2758��ѹ����ڵ����Ҳ����ڶ����ڵ㣬��ô������ѹ��
            CD_node=find(Clique{i}==C(j));
            Clique{i}(CD_node)=[];%���ڶ����ڵ�Ӿֲ������г�ȥ
            C(j)=[];
        else %�����2758��ѹ����ڵ����ܹ��ҵ��ڶ����ڵ㣬��ô������ѹ��
            Matrix(:,vertex_min)=Matrix(:,vertex_min)+Matrix(:,vertex_max);
            Matrix(vertex_min,:)=Matrix(vertex_min,:)+Matrix(vertex_max,:);
            Matrix(vertex_min,vertex_min)=Matrix(vertex_min,vertex_min)-Matrix(vertex_min,vertex_max)-Matrix(vertex_max,vertex_max);
            Matrix(vertex_max,:)=[];
            Matrix(:,vertex_max)=[];
            index=find(D==C(2));
            D(index)=[];
            C(j)=[];
        end
    end
end
end

function f=inherit(ParetoFront1,Clique)

RANK=cell2mat(cellfun(@(S)S(randi(length(S))),Clique,'UniformOutput',false));
% �õ�ÿ���ֲ����ŵ�ģ������
RANK=ParetoFront1(RANK);
% ���ǵ���Щ�ֲ���������ͬһ��ģ�飬��������ѭ������ѡ��ÿ�������ֲ����ŵ�ģ�飬
for j=1:max(RANK)
    % ��ȡ���ڵ� j ��ģ��ľֲ����ŵ�������
    O=find(RANK==j);
    % ��� O �ǿ�
    if length(O)>0
        % �� O �ĵ�һ��Ԫ��ת������ O �����
        O=[O O(1)];
        O(1)=[];
        % ��Ϊ f Ϊ�ھӱ�����Ľ��������������ʹ������ͬһ��ģ��ľֲ����Ż���Ϊ����
        f(find(RANK==j))=O;
    end
end

end