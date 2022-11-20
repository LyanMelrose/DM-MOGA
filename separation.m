function  [A,B,c]= separation(index,Adj_mat)    %c判断是否分裂
M=Adj_mat(index,index);
degree=sum(M,1);
D=find(degree>2);
%% 找到4派系，之后往里面加点

 [ Community,cliques_4 ] = PreCluster(M);



if cliques_4==1
    K=Community;
    [~,A]=sort(degree);
    A=fliplr(A);
%     A= setdiff(A,K); 

visited=zeros(1,length(index));
K1_count=0;
visited(K)=1;
    for i=1:length(A)
       l=degree(A(i));
       l_K=length(find(M(A(i),K)==1));
       if l_K/l>0.5
           K1=[K A(i)];
           visited(A(i))=1;
           K1_count=1;
       end
    end
   if K1_count==1
       K=unique(K1);
   end
  A= index(find(visited==1));
  B=index(find(visited==0));
  C=[];
  if length(B)>0
  if B(1)<A(1)
      C=B;
      B=A;
      A=C;
  end
  end
  if length(A)==length(index)
       c=0; 
  else 
      c=1;
  end
end

  
    
    
    
  
    
    
    






%% 如果没有4派系，找到3派系（从度的小的找），往里面加点
if cliques_4==0
   
D=find(degree==2);
if length(D)==0
    D=find(degree==3);
end
index_i=0;
visited=zeros(1,length(index));
min_degree=inf;
for i=1:length(D)
    K=[D(i) find(M(D(i),:)==1)];
        degree_K=sum(degree(K),2);
       if min_degree>degree_K
         min_degree=degree_K;
         index_i=i;
       end
end
K=[D(index_i) find(M(D(index_i),:)==1)];
visited(K)=1;
[~,A]=sort(degree);
A=fliplr(A);

   for i=1:length(A)
       l=degree(A(i));
       l_K=length(find(M(A(i),K)==1));
       if l_K/l>0.5
           K=[K A(i)];
           visited(A(i))=1;
       end
   end
   K=unique(K);
  A= index(find(visited==1));
  B=index(find(visited==0));
  C=[];
  if length(B)>0
  if B(1)<A(1)
      C=B;
      B=A;
      A=C;
  end
  end
  if length(A)==length(index)
       c=0; 
  else 
      c=1;
  end
end









end


