function QQ=ff(AA,adj_mat,ll,N)
ll=2*ll;
n=length(AA);
nn=length(N);
numVar=size(adj_mat,1);
degree=sum(adj_mat,1);
sum1=0;
index1=[];
for i=1:nn
    if(adj_mat(AA(1),N(i))==1&&adj_mat(AA(2),N(i))==1)
        index1=[index1 N(i)];
    end
end
for j=1:length(index1)
    k=degree(index1(j));
    sum1=sum1+1/log(k);
end
    
QQ=sum1;
% sum2=0;
% for i=1:n
%     for j=1:n
%         
%        k=adj_mat(AA(i),AA(j))-degree(AA(i))*degree(AA(j))/ll;
%         sum2=sum2+k;
%        
%     end
% end
% QQ=sum1+sum2;

% for i=1:n
%     for j=1:n
%         if i~=j&&adj_mat(i,j)==1
%             sum=sum+1;
%         end
%             
%     end
%     
    
    
end
% sum=sum/2;
% m=size(adj_mat,1);
% 
% for i=1:m
%     B(i)=i;
% end
% C=setdiff(B,A);
% sum1=0;
% for kk=1:length(A)
% for k=1:length(C)
%     if(adj_mat(C(k),A(kk))==1)
%         sum1=sum1+1;
%     end
% end
%     
% end

