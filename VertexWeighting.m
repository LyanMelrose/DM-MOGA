function Weights=VertexWeighting(AdjacentMatrix)
%%
numVertex=size(AdjacentMatrix,1);
Weights=zeros(1,numVertex);
for i=1:numVertex
    %% might error
%     neighbors=find(AdjacentMatrix(i,:)==1);
%     neighborsAdjacentMatrix=AdjacentMatrix(neighbors,neighbors);
%     degree=sum(neighborsAdjacentMatrix);
%     maxK=max(degree);
%     CoreVertex=find(degree==maxK);
%     numCore=length(CoreVertex);
%     HighestKCore=neighborsAdjacentMatrix(CoreVertex,CoreVertex);
%     EdgeNum=sum(sum(HighestKCore))/2+numCore;
%     CoreClusterCoefficient=(2*EdgeNum)/((numCore+1)*numCore);
%     Weights(1,i)=maxK*CoreClusterCoefficient;

%% rewrite mcode weighting
    neighbors=find(AdjacentMatrix(i,:)~=0);
    neighborsAdjacentMatrix=AdjacentMatrix([i,neighbors],[i,neighbors]);
%     neighborsAdjacentMatrix=AdjacentMatrix(neighbors,neighbors);
    [HighestKCore,Kmax]=FindHighestKCore(neighborsAdjacentMatrix);
    numCore=size(HighestKCore,2);
    numEdge=sum(sum(HighestKCore))/2;
    CoreClusterCoefficient=2*numEdge/(numCore*(numCore-1));
    Weights(1,i)=Kmax*CoreClusterCoefficient;
end
end

function [HighestKCore,Kmax]=FindHighestKCore(neighborsAdjacentMatrix)
    vertex=1:size(neighborsAdjacentMatrix,2);
    HighestKCore=neighborsAdjacentMatrix(vertex,vertex);
    Kmax=min(sum(HighestKCore));
    while ~isempty(vertex)
        KCore=neighborsAdjacentMatrix(vertex,vertex);
        degree=sum(KCore);
        K=min(degree);
        if K>Kmax
            Kmax=K;
            HighestKCore=KCore;
        end
        DeleVertex= degree==K;
        vertex(DeleVertex)=[];
        if size(vertex,2)-1<=Kmax
            break
        end
    end
end
