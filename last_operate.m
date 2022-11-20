function last_operate(string_m,AdjacentMatrix,e,degree1,Node3)

global  Thirty_Run_minWW
% metric = load('results1\last_operate\MODPSO1_0.35_metrics1.txt');
%  for k=1:10
strNetwork=string_m;
% string_m='0.45';
% e=1;
path = sprintf('results1/%s/ParetoFront/MODPSO1_%s_PF%d.txt',string_m,string_m,e);
ParetoFront=load(path);

% AdjacentMatrix = load('GNExtend\0.45.txt');
% Datalabel=load('GNExtend\real0.45.txt');
% AdjacentMatrix = load('RealWorld\football.txt');
% Datalabel=load('RealWorld\real_label_football.txt');
Adj_mat=AdjacentMatrix;
% m=find(B==1);
% n=find(B==4);
% B(m)=B(m)+3;
% B(n)=B(n)-3;
% V=Datalabel-B;
% index=find(V~=0);
%[Node1,matrix,degree1,edgeslist]=Get_Cliques(Adj_mat,Datalabel);
numVar=size(Adj_mat,1);

V=length(degree1);

for in=1:size(ParetoFront,1)
    A=ParetoFront(in,1:V);
    ParetoFront1(in,1:numVar)= decode1(A,Node3);
    change_node{in}=[];
    %% %%后处理：2种方式-最大的Weighted Wscore
    label= ParetoFront1(in,1:numVar);
    for i=1:length(degree1)
        if(length(Node3{i})>=1)
            for j=1:length(Node3{i})
                
                m=Node3{i}(j);
                if m==114
                    wo=1;
                end
                k=label(m);
                index_last=find(label==k);
%                 if(fff(A,Adj_mat,ll)>fff(index_last,Adj_mat,ll))
                    neighbors=find(Adj_mat(m,:));
                    com_max=Max(neighbors,label);
                    ParetoFront1(in,m)=com_max;
                    if com_max~=k
                        change_node{in}=[change_node{in} m];
                    end
            end
        end    
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/8/16%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [~,index1]=min(ParetoFront(:,V+1));
%     fprintf('minWW :    %g\n',ParetoFront(index1,V+1));
%     fprintf('Paretosize = %g\n',size(ParetoFront1,1));
%     path = sprintf('results1/%s/metrics/MODPSO1_%s_metrics%d.txt',strNetwork,strNetwork,e);
%     savedata1(path,metrics);
    Thirty_Run_minWW=ParetoFront(index1,V+1);
    
end
