function f = value(x,AdjMatrix,ll)

f=[];

%clu_assignment = decode(x);



f(1) = value_1(AdjMatrix,x,ll)/2;
f(2) = value_2(AdjMatrix,x,ll);

end