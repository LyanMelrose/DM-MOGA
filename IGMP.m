

function P=IGMP()
global numVar  AdjacentMatrix
 Adj=AdjacentMatrix;
 popnum=numVar;
 
 

t=5;
N=size(Adj,1);

%A是一个种群
A=[];

%初始化
for i=1:popnum
    A(i)=unidrnd(numVar);
end;

B=A;
   
    for j=1:t
        
        Q1=Modularity(A);
        for k=1:N
            
            neighbors=find(Adj(k,:));
            com_max=Max(neighbors,A);
           B(k)=com_max;
        end
        
        Q2=Modularity(B);
        if Q2>Q1
            A=B;
        else
            break;
        end;
    end

P=A;
end
