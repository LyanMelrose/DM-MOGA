function  C_1_node_j=compute_C_1(j,AdjMatrix,degree)  %%����Լ��ľ���ϵ��

     k=degree(j);
     E=0;
     
     neighbours=find(AdjMatrix(j,:)==1);
            for m=1:k
              for n=1:k
                  if(AdjMatrix(neighbours(m),neighbours(n))==1)
                      E=E+1;
                  end
              end
            end
            if k>1
            C_1_node_j=E/(k*(k-1));
            else
                C_1_node_j=1;
            end
end