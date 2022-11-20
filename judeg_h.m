function pos=judeg_h(AA,BB)
A=AA;
B=BB;
 n=length(A);
 m=length(B);
 sum1=0;
 for i=1:n
     for j=1:m
         if(A(i)==B(j))
             sum1=sum1+1;
         end
     end
 end
 
k=min(n,m);
if(sum1/k>=0.5)
    pos=1;
else 
    pos=0;
end