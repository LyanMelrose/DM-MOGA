function Node=GET_CK()
global AdjMatrix degree
%AdjMatrix=load('GNExtend/0.4.txt');
AdjMatrix=[0 1 0 0 0 0 0;
           1 0 1 0 1 0 1;
           0 1 0 1 1 1 0;
           0 0 1 0 1 0 1;
           0 1 1 1 0 1 1;
           0 0 1 0 1 0 1;
           0 1 0 1 1 1 0;];
           
degree=sum(AdjMatrix,1);
%Datalabel=load('GNExtend/real0.4.txt');
numVar=size(AdjMatrix,1);

for i=1:numVar
    N{i}=find(AdjMatrix(i,:)~=0);     %%是个变动的数组
    S{i}=i;
end
visited=zeros(1,numVar);



%外面还要再加层循环 i：表示团号
    while(1)
           for j=1:length(N{i})
               node_j=N{i}(j);             
               Over_half=Over(S{i},node_j);
               if Over_half>0.5                    %%大于一半的加入
                   S{i}=[S{i} node_j];
                    N{i}=unique([N{i} N{j}]); 
                    
                     N{i}=setdiff(N{i},S{i});
               end
           end
         
           C_i=compute_Ci(S{i});
           
           if C_i==1
               s{i}=[S{i} N{i}];
               N{i}=unique([N{i} N{j}]);           %%构成一个派系并入
               continue;
           end
           
           for j=1:length(N{i})
               node_j=N{i}(j);
               C_1_node_j=compute_Ci(node_j);
               if C_1_node_j==1
                   S=[S{i} N{j}];                  %%属于一个派系的一点
                   N{i}=unique([N{i} N{j}]); 
                   break;
               end
           end
           if j~=length(N{i})
               continue;
           end
           
       
           C_i=compute_Ci(S{i});
           for j=1:length(N{i})
               node_j=N{i}(j);
               C_1_node_j=compute_C_1(node_j);
               C_i=compute_Ci(S{i});

               if C_i>C_1_node_j
                   del_C(j)=-inf;
               else
                   C_2_node_j=compute_C_2(S{i},node_j);   %%%函数没写
                   del_C(j)=C_1_node_j- C_2_node_j;                   %%%下次接着写
               end
           end
           
           count=0;
           max_del_C=-inf;
           pos=0;
           for j=1:length(N{i})
               node_j=N{i}(j);
               if del_C(j)<0
                   count=count+1;
               else del_C(j)>=max_del_C
                   max_del_C=del_C(j);
                   pos=j;
               end
           end
           
           if count<length(N{i})
               S{i}=[S{i} N{i}(pos)];
               
               N{i}=unique([N{i} N{N{i}(pos)}]);
               N{i}=setdiff(N{i},S{i});
               continue;
           else length(N{i})==1
              continue;
           end
           
       end

end












