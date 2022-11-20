function gene1 = IGLP(numVar,node,Node3)

%

% test
% X = load('karateAdjMatrix');
% AdjacentMatrix=X.karateAdjMatrix;% karate为一个结构体，不能直接赋值给AdjMatrix
% gene=[1 1 1 1 1 1 1 1 2 2 1 1 1 1 2 2 1 1 2 1 2 1 2 2 2 2 2 2 2 2 2 2 2 2];
% numVariables=34;
% for i=1:numVariables
%     node(i).neighbours=find(AdjacentMatrix(i,:)==1);
%     node(i).degree=size(node(i).neighbours,2);
% end


gene = 1:numVar;


for n = 1 : 5  
    a=randperm(numVar);
    for  i = 1 : numVar    
        
        NeighborSize = length(node(a(i)).e);
        if NeighborSize == 0   
            %			gene(i) = 0;            
        else            
            if NeighborSize == 1                
                gene(a(i)) = gene(node(a(i)).e(1));                
            else                
%                 sum = 0;
                maxr = -1;%//record index of i's neighbour which ...
                label = -1;
                temp = 1;                
                for j = 1 : NeighborSize                   
                    counter = 1; %//record no. of nodes that has same label with j 
                    for k = j + 1 : NeighborSize 
                        p = gene(node(a(i)).e(j));
                        q = gene(node(a(i)).e(k));
                        if p == q
                            counter = counter + 1;
                        end
                    end %//end k
                    if temp < counter                       
                        maxr = j;
                        temp = counter;
                    end
                end %//end j
                for l =1 : NeighborSize                    
                    u = gene(node(a(i)).e(l));
                    v = gene(i);
                    if  u == v                        
                        label = u;%%找到与当前节点相邻的点在同一个团，不是孤立点
                    end
                end %//end l
                if label ~= -1 && maxr == -1                    
                    gene(a(i)) = label;                   
                else                    
                    if maxr ~= -1                        
                        gene(a(i)) = gene(node(a(i)).e(maxr)); 
                        if(length(Node3{a(i)})>=2) gene(a(i))=a(i); end
                    else                       
                        randneighbor = randi(NeighborSize);
                        %randneighbor = 2;%test
                        gene(a(i)) = gene(node(a(i)).e(randneighbor));  
                        if(length(Node3{a(i)})>=3) gene(a(i))=a(i); end
                    end
                end
            end % if NeighborSize == 1 
        end % if NeighborSize == 0
    end %//end i
end %//end n
gene1=gene;
% for i=1:numVar
% if(length(node(i).e)>=3)
%     gene1(i)=i;
% end

% gene=decode1(gene,Node3);
% 
% Q=Modularity(gene);
% if Q~=0
%     gene(numVar+1)=Q;
%     break;
% end
end