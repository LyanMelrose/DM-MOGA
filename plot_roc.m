function AUC = plot_roc( predict, ground_truth )
%  predict       - �������Բ��Լ��ķ�����
%  ground_truth - ���Լ�����ȷ��ǩ,����ֻ���Ƕ����࣬��0��1
%  auc            - ����ROC���ߵ������µ����

% �����ground_truth������������Ŀpos_num�͸���������Ŀneg_num
pos_num = sum(ground_truth==1);
neg_num = sum(ground_truth==0);

% ���Լ��������ĸ���
m=size(ground_truth,1);
% ��������
[pre,Index]=sort(predict);
% ��ground_truth��������ΪԪ����������predictһһ��Ӧ��
ground_truth=ground_truth(Index);
% ��ʼ��x��y
x=zeros(m+1,1);
y=zeros(m+1,1);
AUC=0;
% ����ʼ����Ϊ��1.0, 1.0��
x(1)=1;y(1)=1;

% ����ȡ��һ������������ROC���ߺ�AUCֵ��
for i=2:m
TP=sum(ground_truth(i:m)==1);FP=sum(ground_truth(i:m)==0);
x(i)=FP/neg_num;
y(i)=TP/pos_num;
AUC=AUC+(y(i)+y(i-1))*(x(i-1)-x(i))/2;
end
% ���յ㱻��Ϊ(0,0)
x(m+1)=0;y(m+1)=0;
% ����AUCֵ
AUC=AUC+y(m)*x(m)/2;

% ��ͼ
plot(x,y,'lineWidth',1.5);
xlabel('False Positive Rate');
ylabel('True Positive Rate');
end
