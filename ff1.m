function QQ=ff1(AA,adj_mat,degree)
%adj_mat=load('RealWorld\karate.txt');
%AA=[ 1 12];
BB=adj_mat(AA,AA);
degree_BB=sum(BB,1)
k_in=sum(degree_BB)/2;
k_out=sum(degree(AA))-k_in*2;
QQ=k_in/(k_in+k_out);
end


% n=length(AA);
% numVar=size(adj_mat,1);
% degree=sum(adj_mat,1);
% sum1=0;
% index1=AA;
% % for i=1:numVar
% %     if(adj_mat(AA(1),i)==1&&adj_mat(AA(2),i)==1)
% %         index1=[index1 i];
% %     end
% % end
% % for j=1:length(index1)
% %     k=1/degree(index1(j));
% %     sum1=sum1+k;
% % end
% %     
% % QQ=sum1
% sum2=0;
% for i=1:n
%     for j=1:n
%         
%        k=adj_mat(AA(i),AA(j))-degree(AA(i))*degree(AA(j))/ll;  %%偏度数低的方向走
%        
%         sum2=sum2+k;
%        
%     end
% end
% QQ=sum2;
% 
% % for i=1:n
% %     for j=1:n
% %         if i~=j&&adj_mat(i,j)==1
% %             sum=sum+1;
% %         end
% %             
% %     end
% %     
%     
%     
% end
% % sum=sum/2;
% % m=size(adj_mat,1);
% % 
% % for i=1:m
% %     B(i)=i;
% % end
% % C=setdiff(B,A);
% % sum1=0;
% % for kk=1:length(A)
% % for k=1:length(C)
% %     if(adj_mat(C(k),A(kk))==1)
% %         sum1=sum1+1;
% %     end
% % end
% %     
% % end
% 
