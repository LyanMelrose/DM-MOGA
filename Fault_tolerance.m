function   [ Node3,chromosomes, degree1, edgeslist,V,fault_node]= Fault_tolerance(Node3,chromosomes, degree1, edgeslist,V,fault_node,AdjMatrix)
popsize=size(chromosomes,1);
numVar=size(AdjMatrix,1);
BB=chromosomes(:,[V+1 V+2]);   %找出parate面上第一个曲面上的个体
B=P_sort(BB,'all');
B=find(B==1);

BB=[];                        %给parate第一层上的解译码
for l=1:length(B)
BBB=decode(chromosomes(B(l),1:V));
BB=[BB; BBB];
end

g=ones(1,numVar);             %找出错误相邻的点，用来查找散点和散点们处于哪个团，之后看fault点所在的团的标签和散点标签做个统计，如果相异的多，那么说明这个团对这个点吸引不大。只是和这个fault
gg=ones(1,numVar);                                                                                                                             %点相连接紧密。
for i=1:V
if(length(Node3(i).e)>=3)
g(Node3(i).e)=0;
end
gg(Node3(i).e)=i;
end


remove=[];
for i=1:length(fault_node)
    A=find(AdjMatrix(fault_node(i),:)==1);
    lable=BB(gg(fault_node(i)));
    Matrix=BB(:,gg( A(find(g(A)==1))));
    [l,r]=size(Matrix);
    count=0;
    for j=1:l
        for jj=1:r
            if Matrix(j,jj)~=lable
                count=count+1;
            end
        end
    end
    if count/(l*r)>0.25
        remove=[remove fault_node(i)];
    end
end



%%%更新 Node3,chromosomes, degree1, edgeslist,V。。。。。。。。。。
for i=1:length(remove)
    a=gg(remove(i));
    Node3(a).e=setdiff(Node3(a).e,remove(i));
    Node3(V+i).e=remove(i);
end
t=length(remove)+V;
lll=1;
for m=1:t                              %%社团第一个标号排序，以利于matrix折叠
    B1(m)=Node3(m).e(1);
end
[~,C1]=sort(B1);
for lll=1:t
    Node3(lll).e=Node3(C1(lll)).e;
end
matrix=AdjMatrix;
numVar=size(AdjMatrix,1);
times=0;
for z=1:numVar
    BBB(z)=z;
end
D=BBB;
for i=1:t
    C=Node3(i).e;
     vertex_min=find(D==C(1));
  while length(C)>1
       j=2;
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
%%%测试折叠的正确性
for mm=1:t
   Node3(mm).e
   Datalabel(Node3(mm).e)
end
edgeslist = edges_list(matrix,t,Node3);
degree=sum(AdjMatrix,1);
for i=1:t
   degree1(i)=sum(degree(Node3(i).e),2);
end

%%对种群进行更新
lable=zeros(1,t+2);
chromosomes_1=zeors(popsize,t+2);
d=1;
for i=1:t
    if(Node3(i).e(1)==remove(d))
        d=d+1;
        lable(i)=0;
    else
        lable(i)=i-d+1;
    end
end
% for i=1:popsize
%     for j=1:t
%         if lable(j)==0
%             R=randperm(edges_list(j).n);
%             chromosomes_1(i,j)=edges_list(j).e(R(1));
%         else
%         chromosomes_1(i,j)=chromosomes(i,lable(j));
%         end
%         
%         
%     end
%       f(i,V + 1: K) = evaluate_objective(f(i,1:V),ll,degree1);
%       chromosomes_1(i,t+1:t+2)=evaluate_objective(chromosomes_1(i,1:t),ll,degree1)
% end
% idealp = min(chromosomes_1(:,t+1:t+2));






end








    
    
    
        
    
    
    
    
    
    
    

    

















