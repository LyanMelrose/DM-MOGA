function Over_half=Over(A,i)
global AdjMatrix degree
count_V=0;
for j=1:length(A)
    if(AdjMatrix(A(j),i)==1)
        count_V=count_V+1;
    end
end
Over_half=count_V/degree(i);
end