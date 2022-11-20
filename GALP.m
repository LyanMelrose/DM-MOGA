function chorm=IGLP()
global node   numVar
chorm=zeros(1,numVar);
for ii=1:numVar
    chorm(ii)=ii;
end
for i=1:numVar
    index=randperm(node(i).degree);
    x=index(1);
    chorm(i)= chorm(node(i).neighbours(x));

end