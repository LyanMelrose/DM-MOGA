function [Node3,matrix,degree1,edgeslist]=Get_Cliques(Adj_mat,Datalabel)
global Node3 edge_num numVar 
numVar=size(Adj_mat,1);
adj_mat=Adj_mat;
[V,V] = size(adj_mat); 
degree=sum(adj_mat,1);
edge_num=sum(degree)/2;
ll=edge_num;
%% 计算点的邻居集和共同邻居 ——node(i)表示i点的邻居集——
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
%% 点合并成团--采用相似度最大的2个点合并
for i=1:numVar
    if i==57
        stop=1;
    end
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
Node(i).e=[Node(i).e index];
      end
    A=[];
    end
end
%% 对numVar个团进行合并--t为合并之后的团数目
visited=zeros(1,numVar);
for m=1:numVar
    A=Node(m).e;
   
     if visited(1,m)==0
         while(1)
             Dep=length(Node(m).e);
     for n=1:numVar
         
         if m~=n&&visited(1,n)==0
           B=Node(n).e;
          
           k=common(Node(m).e,B);
           
          if k==1
              if m==4
                  cc=0;
              end
            Node(m).e=unique([Node(m).e B]);
            visited(1,n)=1;
              
          end
         end
     end
     Dep_1=length(Node(m).e);
     if Dep_1==Dep
         break;
     end
        end
     end
end
 t=1;
for kk=1:numVar
   if visited(1,kk)==0
       if t==13
           kk;
       end
        Node1(t).e=Node(kk).e;
        Node1(t).e;
        Datalabel(Node1(t).e);
        t=t+1;
    end
end
t=t-1;
%% 测试与真实情况的对比和合并有无出错的情况+计算初步划分得到的Q值---计算得到的Q<0.4的都需要拆分
for o=1:t
   o;
    Node1(o).e;
    Datalabel(Node1(o).e);
end
b=0;
for o=1:t
    b=b+length(Node1(o).e);
end
clu_assignment1=zeros(1,numVar);
for i=1:t
    
        clu_assignment1(Node1(i).e)=i;
end
 Q = compute_dense(Adj_mat,clu_assignment1,edge_num);
 Q_index=find(Q~=0);
 Q=Q(Q_index);
Clique_Q=min(Q);



%% 测试合并的结果对不对+人为的调控错误的团+如合并的效果不好则分裂；  具体怎么拆分？
counter_t=t;
for mm=1:t
%     mm;
%     Node1(mm).e;
%     Datalabel(Node1(mm).e);
    if length(Node1(mm).e)>4
    [A,B,c]= separation(Node1(mm).e,Adj_mat);
  
    if c==1
        Node1(mm).e=A;
        counter_t=counter_t+1;
        Node1(counter_t).e=B;
    end
    end
    
end
a=[];
t=counter_t;
for o=1:t
   o;
   if length(find(Node1(o).e==981))==1;
       a=[a o];
   end
   
    Datalabel(Node1(o).e)
    
end

for m=1:t
    B1(m)=Node1(m).e(1);
end
[~,C1]=sort(B1);
for lll=1:t
    Node3(lll).e=Node1(C1(lll)).e;
end
% Node1(25).e=Node1(21).e([5 6 7 8 9 10])
% temp=Node1(21).e([1 2 3 4 11]);
% Node1(21).e=[];
% Node1(21).e=temp;
% Node1(26).e=Node1(14).e([3 6 7 8]);
% temp1=Node1(14).e([1 2 4 5]);
% Node1(14).e=[];
% Node1(14).e=temp1;
% 
% Node1(27).e=Node1(13).e([3 4 5 7 9]);
% temp2=Node1(13).e([1 2 6 8 10]);
% Node1(13).e=[];
% Node1(13).e=temp2;
% 
% Node1(28).e=Node1(16).e([2 3 5 6 8]);
% temp3=Node1(16).e([1 4 7]);
% Node1(16).e=[];
% Node1(16).e=temp3;
%% 得到合并的团所对应的邻接矩阵--matrix---D表示要合并过点，合并完的点从D里面去掉
matrix=Adj_mat;
times=0;
for z=1:numVar
    BBB(z)=z;
end
D=BBB;
for i=1:t
    if i==86
        lllll=0;
    end
    C=Node3(i).e;
     vertex_min=find(D==C(1));
  while length(C)>1
       j=2;
       if C(2)==951

           
           jj=0;
       end
      
       vertex_max=find(D==C(2));
       if length(vertex_max)==0
           CD_node=find(Node3(i).e==C(j));
           Node3(i).e(CD_node)=[];
           C(j)=[];
       else
           
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
end
%% 团的标号必须和matrix对应+测试部分
% t=length(matrix);
% for m=1:t
%     B1(m)=Node1(m).e(1);
% end
% [~,C1]=sort(B1);
% for lll=1:t
%     Node3(lll).e=Node1(C1(lll)).e;
% end
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

for mm=1:t
   Node3(mm).e;
   Datalabel(Node3(mm).e);
end
%% 得到信息--edgeslist表示团的邻接团集--deree1表示团的度
%%有重叠点，调试下
edgeslist = edges_list(matrix,t,Node3);
degree=sum(adj_mat,1);
for i=1:t
degree1(i)=sum(degree(Node3(i).e),2);
end
end