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