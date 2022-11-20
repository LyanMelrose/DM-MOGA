function [Community,MaxDeltaQ] = FastCommunity1()
adj_mat = load('RealWorld\125_point.txt');
Datalabel=load('RealWorld\real_label_125_point.txt');
numVar=size(adj_mat,1);
[V,V] = size(adj_mat); 

Adj_mat=adj_mat;

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
  
           for j=1:length(A)
           B=[i A(j)];
           del_Q=ff(B,adj_mat,ll);
       
            if del_Q>maxQ
               index=A(j);
               maxQ=del_Q;
            end
           end
%     [~,index]=max(Q(i).e);
%     A(index)=[];
%     Q(i).e(index)=[];
%     [~,index]=max(Q(i).e);
%     second_max=A(index)
    Node(i).e=[Node(i).e index];
        end
    A=[];
    end
end
%%ºÏ²¢
for mmm=1:numVar
    Node(mmm).e
end
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

        


for mm=1:(t-1)
    mm
    Datalabel(Node1(mm).e)
end

   
end
function logicall=common(A,B)
logicall=0;
[~,a]=size(A);
[~,b]=size(B);
for i=1:a
    for j=1:b
        if A(i)==B(j)
            logicall=1;
            break;
        end
    end
end

end
























