function [Indicators,Avg_AUC] = FiveFold_SVM(data,label)
%Ҫ�����ݼ�ÿ�д���һ��������labelΪһ�д洢��������ǩ
Kfold=5;
% M:���������� N:һ��������Ԫ������
[M,N]=size(data);
% ��������������ְ���indicesΪ����M��Ԫ�ص�һά�������ڶ�������Ϊlabelû�����⡣
% ÿ��Ԫ�ض�Ӧһ�����������д洢��[1,5]��������Ϊ�������ķ����š�
indices=crossvalind('Kfold',label,Kfold);

%������֤k=5��5����������Ϊ���Լ�
AUC=zeros(Kfold,1);
for k=1:Kfold
    % ���test��Ԫ�������ݼ��ж�Ӧ�ĵ�Ԫ���
    test = (indices == k);
    % train��Ԫ�صı��Ϊ��testԪ�صı��
    train = ~test;
    
    % ͬ���أ�������������Ҫ������ľ��󡪡���Ϊ��������Ϊ������
    model = fitcsvm(data(train,:),label(train,:));
    test_predict = predict(model,data(test,:));
    %������ȷ��
    Indicators(k,1) = EvaluationCriteria_KnownClass(label(test),test_predict);
    AUC(k,1) = plot_roc(test_predict,label(test));
end
accuracy=mean([Indicators(:,1).accuracy]);
recall=mean(mean([Indicators(:,1).recall]));
precision=mean(mean([Indicators(:,1).precision]));
F_measure=mean(mean([Indicators(:,1).F_measure]));
G_mean=mean([Indicators(:,1).G_mean]);
Indicators=struct('accuracy',accuracy,'recall',recall,'precision',precision,'F_measure',F_measure,'G_mean',G_mean);
Avg_AUC=mean(AUC);
end
