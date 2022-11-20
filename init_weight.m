function [weights neighbors] = init_weight(popsize, niche)
% init_weights and neighbors.

%ÿ����������Ȩ�أ��ֱ��ʼ��Ϊ������/Ⱥ���������Ⱥ�����-i��/Ⱥ�������
weights = [];
for i=0:popsize-1
    weight=zeros(1,2);
    weight(1)=i/(popsize-1);
    weight(2)=(popsize-i-1)/(popsize-1);
    weights = [weights;weight];
end


%Set up the neighbourhood.��ʼ��ÿ��������ھӡ�
leng=size(weights,1);
distanceMatrix=zeros(leng, leng);
for i=1:leng
    for j=i+1:leng
        % �������и����Ȩ�أ��õ������ľ������
        A=weights(i,:)';B=weights(j,:)';
        distanceMatrix(i,j)=(A-B)'*(A-B);
        distanceMatrix(j,i)=distanceMatrix(i,j);
    end
    % ����ÿ�����壬���ݸø�������������ľ�����������������������
    [s,sindex]=sort(distanceMatrix(i,:));
    % ѡȡ��С��niche��������Ϊ�ø�����ھӡ�
    neighbors(i,:)=sindex(1:niche);
end
   
end