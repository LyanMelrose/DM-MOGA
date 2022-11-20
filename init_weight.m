function [weights neighbors] = init_weight(popsize, niche)
% init_weights and neighbors.

%每个个体有俩权重，分别初始化为个体标号/群体个数、（群体个数-i）/群体个数。
weights = [];
for i=0:popsize-1
    weight=zeros(1,2);
    weight(1)=i/(popsize-1);
    weight(2)=(popsize-i-1)/(popsize-1);
    weights = [weights;weight];
end


%Set up the neighbourhood.初始化每个个体的邻居。
leng=size(weights,1);
distanceMatrix=zeros(leng, leng);
for i=1:leng
    for j=i+1:leng
        % 根据所有个体的权重，得到个体间的距离矩阵。
        A=weights(i,:)';B=weights(j,:)';
        distanceMatrix(i,j)=(A-B)'*(A-B);
        distanceMatrix(j,i)=distanceMatrix(i,j);
    end
    % 对于每个个体，根据该个体与其他个体的距离对其他个体进行升序排序。
    [s,sindex]=sort(distanceMatrix(i,:));
    % 选取最小的niche个个体作为该个体的邻居。
    neighbors(i,:)=sindex(1:niche);
end
   
end