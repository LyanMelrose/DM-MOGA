function  C_2_node_j=compute_C_2(A,j,AdjMatrix,degree)

    
     E=0;
     neighbours=find(AdjMatrix(j,:)==1);
   
     B=setdiff(neighbours,A);
     k=length(B);
     
            for m=1:k
              for n=1:k
                  if(AdjMatrix(B(m),B(n))==1)
                      E=E+1;
                  end
              end
            end
          
            if k>1
                C_2_node_j=E/(k*(k-1)); 
            elseif k==1
                C_2_node_j=1;
            elseif k==0
                C_2_node_j=0;
            end
            
                
                


end