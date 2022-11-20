function  C_i=compute_Ci(A,AdjMatrix,degree)
   
        sum_C=0;
        E=0;
        for i=1:length(A)
            j=A(i);
            k=degree(j);
            neighbours=find(AdjMatrix(j,:)==1);
            for m=1:k
              for n=1:k
                  if(AdjMatrix(neighbours(m),neighbours(n))==1)
                      E=E+1;
                  end
              end
            end
            if k==1
                m=1;
            else
                m=E/(k*(k-1));
            end
            sum_C=sum_C+m;   %--------当k=1时候怎么处理这个公式的。
            E=0;   
        end
            C_i=sum_C/length(A);
   end
