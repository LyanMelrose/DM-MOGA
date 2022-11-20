clc;
clear;
tic;

global numObjectives popsize numVar  max_gen mutate_posibility M niche adj_mat V ll pc pg pm
global   strNetwork Datalabel idealp weights neighbors  edgeslist matrix degree
global Thirty_Run_maxQ Thirty_Run_maxNMI   Adj_mat Node AdjacentMatrix Node3

max_gen = 300;
M = 2;
pg=1;
pm=0.3;
pc=1;
runtimes=10;
strNetwork='0.35';
AdjacentMatrix = load('Net_31773_HardThreshold.mat');
% Datalabel=load('GNExtend\real0.35.txt');


Adj_mat=int8(AdjacentMatrix);
clear global AdjacentMatrix;




% Adj_mat=[0 1 0 1 0 0 0 0 0 0 0 0 0 0 0;1 0 1 1 0 0 0 0 0 0 0 0 0 1 0;0 1 0 1 0 0 0 0 0 1 0 0 0 0 0;
% 1 1 1 0 0 0 0 0 1 0 0 0 0 0 0;0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 ;0 0 0 0 0 0 1 1 1 0 0 0 0 0 0;
% 0 0 0 0 1 1 0 0 1 0 0 0 0 0 0;0 0 0 0 1 1 0 0 1 0 1 0 0 0 0;0 0 0 1 1 1 1 1 0 0 0 0 0 0 0;
% 0 0 1 0 0 0 0 0 0 0 1 0 1 1 1;0 0 0 0 0 0 0 1 0 1 0 1 0 1 1;0 0 0 0 0 0 0 0 0 0 1 0 1 1 1;
% 0 0 0 0 0 0 0 0 0 1 0 1 0 1 1;0 1 0 0 0 0 0 0 0 1 1 1 1 0 0;0 0 0 0 0 0 0 0 0 1 1 1 1 0 0];

% Adj_mat=[0 1 1 1 0 0 0 0 0;
%          1 0 1 0 0 0 0 0 0;
%          1 1 0 1 0 0 0 0 0;
%          1 0 1 0 1 1 0 0 0;
%          0 0 0 1 0 1 1 1 0;
%          0 0 0 1 1 0 1 1 0;
%          0 0 0 0 1 1 0 1 1;
%          0 0 0 0 1 1 1 0 0;
%          0 0 0 0 0 0 1 0 0;];


numVar=size(Adj_mat,1);
adj_mat=Adj_mat;
[V,~] = size(adj_mat); 

degree=sum(adj_mat,1);

 ll=0;
for i=1:V
    for j=i:V
        ll=adj_mat(i,j)+ll;
    end
end
for i=1:numVar
     common_neighbours=[];
    node(i).neighbours=find(Adj_mat(i,:)==1);
    
    for k=1:length(node(i).neighbours)
   
        temp=0;
     common_neighbours(k)=0;
    for j=1:numVar
      
        if(Adj_mat(i,j)==1&&Adj_mat(j,node(i).neighbours(k))==1)
         common_neighbours(k)=common_neighbours(k)+1;
         temp=temp+1;
        end
    end
   
        
        
        
    end
    
    node(i).neighbours=[node(i).neighbours;common_neighbours];
 
    node(i).degree=length(node(i).neighbours);
    
end

for i=1:numVar
   A=node(i).neighbours(2,:);
   a=sum(A);
   if a==0
       for j=1:length(A)
           A(j)=0;
       end
   else
   for j=1:length(A)
       A(j)=A(j)/a;
   end
   node(i).neighbours(2,:)=A;
   
    
   end
end
for i=1:numVar
     Node(i).e=i;
    if length(node(i).neighbours(1,:))==1
        Node(i).e=[i node(i).neighbours(1,:)];
   
    else if sum(node(i).neighbours(2,:))==0
            Node(i).e=[i];
        else
           A=node(i).neighbours(1,:);
           maxQ=-inf;
          Q=[];
           for j=1:length(A)
           B=[i A(j)];
           del_Q=ff(B,adj_mat,ll);
           Q=[Q del_Q];
            if del_Q>maxQ
               index=A(j);
               maxQ=del_Q;
            end
           end
% for i=1:numVar
%     if i==10
%         uu=0;
%     end
%     index=0;
%     Node(i).e=i;
%     A=node(i).neighbours(1,:);
%     maxQ=-inf;
%     Q(i).e=[];
%     for j=1:length(A)
%         B=[i A(j)];
%         del_Q=ff1(B,adj_mat,ll);
%         Q(i).e=[Q(i).e del_Q];
%         if del_Q>maxQ
%             index=A(j);
%             maxQ=del_Q;
%         end
%     end
%     [~,index]=max(Q(i).e);
%     A(index)=[];
%     Q(i).e(index)=[];
%     [~,index]=max(Q(i).e);
%     second_max=A(index)
% indexxx=max(Q);
% N=find(Q==indexxx);
%  if length(N)==1
    Node(i).e=[Node(i).e index];
%  end
        end
    A=[];
    end
end
%%合并
% for mmm=1:numVar
%     Node(mmm).e;
% end
visited=zeros(1,numVar);
for m=1:numVar
    A=Node(m).e;
     if visited(1,m)==0
     for n=1:numVar
         if m~=n
        
         if n==11
             o=0;
         end
        B=Node(n).e;
        k=common(Node(m).e,B);
        if k==1
            Node(m).e=unique([Node(m).e B]);
            visited(1,n)=1;
        end
         end
     end
     end
end
 t=1;
for kk=1:numVar
   
    if visited(1,kk)==0
        Node1(t).e=Node(kk).e;
        t=t+1;
    end
end

        


% for mm=1:(t-1)
%     mm;
%     Datalabel(Node1(mm).e);
% end


 


%%以上为对网络初步处理 
%% %%以下为多目标优化团
matrix=Adj_mat;
 %%团所对应的矩阵 
 times=0;
 
 
for z=1:numVar
    B(z)=z;
end
D=B;
m=t-1;
for i=1:m
    C=Node1(i).e;
    index44=find(D<=C(1));
    times_1=size(index44,2);
       vertex_min=times_1;
      
  
 
   while length(C)>1
       j=2;
       index3=find(D<=C(2));
       times=size(index3,2);
        vertex_max=times;
       
    n=size(matrix,1);
    for step=1:n
        matrix(step,vertex_min)=matrix(step,vertex_min)+matrix(step,vertex_max);
        matrix(vertex_min,step)=matrix(vertex_min,step)+matrix(vertex_max,step);
    end;
    matrix(vertex_min,vertex_min)=matrix(vertex_min,vertex_min)-matrix(vertex_min,vertex_max)-matrix(vertex_max,vertex_max);%eii更新必须只要加上eij即可，即加上社团内部的边。
    matrix(vertex_max,:)=[];
    matrix(:,vertex_max)=[]; 
     index4=find(D==C(2));
     D(index4)=[];
    C(j)=[];
   

   end
end
%% 社团初始化
% for i=1:14
% community{i}=Node1(i).e;
% end
% adj2pajek(adj_mat,community,[1],'dolphin','C:\Documents and Settings\Administrator\桌面\网络画图')
t=length(matrix);
for m=1:t
    B1(m)=Node1(m).e(1);
    
end
[~,C1]=sort(B1);
for lll=1:t
    Node3(lll).e=Node1(C1(lll)).e;
end
%%测试团和matrix矩阵是映射关系
% matrix;
% u=0;
% for z=1:length(Node3(4).e)
%     for b=z+1:length(Node3(4).e)
%         if(Adj_mat(Node3(4).e(z),Node3(4).e(b))==1)
%             u=u+1;
%         end
%     end
% end
        

        
edgeslist = edges_list(matrix,m);
degree=sum(adj_mat,1);

for i=1:t
degree1(i)=sum(degree(Node3(i).e),2);
end

        
        
    
    


    
root = sprintf('results1/%s/ParetoFront',strNetwork);  %%建个文件保存实验数据
if ~isdir(root) %判断路径是否存在
    mkdir(root); 
end
root = sprintf('results1/%s/metrics',strNetwork);
if ~isdir(root) %判断路径是否存在
    mkdir(root);
end

Thirty_Run_maxQ=[];
Thirty_Run_maxNMI=[];
for o=1:t
    Node3(o).e
    Datalabel(Node3(o).e)
    
end

for mg=1:runtimes
    fprintf('第%s次运行\n',num2str(mg));

    TMOEAD(matrix,Node1,mg,degree1);

    fprintf('\n');
end 

[~,index1]=max(Thirty_Run_maxQ(:,1));%找到最大的Q索引
[~,index2]=max(Thirty_Run_maxNMI(:,2));%找到最大的NMI索引

fprintf('Thirtyrun_results:\n');
fprintf('maxQ =    %g %g\n',Thirty_Run_maxQ(index1,1),Thirty_Run_maxQ(index1,2));    %输出最大的Q时的Q和NMI
fprintf('maxNMI =  %g %g\n',Thirty_Run_maxNMI(index2,1),Thirty_Run_maxNMI(index2,2));%输出最大的NMI时的NMI和Q

path = sprintf('results1/%s/metrics/MODPSO1_%s_Thirty_Run_maxQ.txt',strNetwork,strNetwork);%保存实验数据取平均值

savedata1(path,[Thirty_Run_maxQ;0 0;mean(Thirty_Run_maxQ(:,1)) mean(Thirty_Run_maxQ(:,2))]); 
path = sprintf('results1/%s/metrics/MODPSO1_%s_Thirty_Run_maxNMI.txt',strNetwork,strNetwork);
savedata1(path,[Thirty_Run_maxNMI;0 0;mean(Thirty_Run_maxNMI(:,1)) mean(Thirty_Run_maxNMI(:,2))]); 
