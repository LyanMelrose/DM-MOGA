function QQ=fff(AA,adj_mat,ll)
sum2=0;
A=adj_mat(AA,AA);
degree=sum(adj_mat,1);
m=0;
n=size(AA,2);
for i=1:n
    for j=1:n
        if adj_mat(AA(i),AA(j))==1
            m=m+1;
        end
    end
end

for i=1:n
    for j=1:n
        
       k=adj_mat(AA(i),AA(j))-degree(AA(i))*degree(AA(j))/2*m;
        sum2=sum2+k;
       
    end
end
QQ=sum2/n;
end


