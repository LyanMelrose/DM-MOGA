function Q = Modularity2222(solution)
%% ��ȡ������Ȩ�����ĳ�ֻ��ֵ�ģ���
%% ģ��ȵ�����ͼ��㷽������������:
% ��Finding and evaluating community structure in network��(2004)
% ����
% AdjacentMatrix - ������ڽӾ���
% clusters - ��������Ż��֣���1���ǽڵ��ţ���2�������ű�ŵı任����3�������ű��
% ���
% Q - ���ֵ�ģ���
% e - k*k����k��������Ŀ������Ԫ�صĶ������������
% a - �������������
% out_clusters - �����뺬����ͬ
global numVar edge_num AdjMatrix;
clusters=zeros(3,numVar);
clusters(1,:)=1:numVar;
clusters(3,:)=solution;

%����Qֵ��Ҫ��������clusters�����е���Ϣ������Ĵ��뽫�������ű����Ϣת���ɾ���ı����Ϣ��ת����Ϣ���ڵڶ���
uni = unique(clusters(3,:));
for i = 1:length(uni)
    idices = find(clusters(3,:) == uni(i));
    clusters(2,idices) = i;
end
degree = sum(AdjMatrix,2); %��AdjMatrix��ÿ����ͣ�degree���ÿ���ڵ�Ķ�
edge_num=sum(degree)/2;
m = edge_num; % ����ıߵ���Ŀ
Q = 0 ;
k = numel(unique(clusters(2,:))); % ������Ŀ
e = zeros(k);
for i=1:k
    idx = find(clusters(2,:)==i);
    labelsi = clusters(1,idx);
    for j=1:k
        idx = find(clusters(2,:)==j);
        labelsj = clusters(1,idx);
        for ii=1:length(labelsi)
            vi = labelsi(ii);
            for jj=1:length(labelsj)
                vj = labelsj(jj);
                e(i,j) = e(i,j)+AdjMatrix(vi,vj);
            end
        end        
    end
end

for i=1:k
    Q = Q + e(i,i)-sum(e(i,:))^2/(2*m);
end
Q=Q/(2*m);

% e = e/2;
% % for i=1:k
% %     e(i,i) = e(i,i)*2;
% % end
% 
% e = e/m;
% % a = [];
% for i=1:k
%     ai = sum(e(i,:));
% %     a = [a; ai];
%     Q = Q + e(i,i)-ai^2;
% end

% out_clusters = clusters;

end