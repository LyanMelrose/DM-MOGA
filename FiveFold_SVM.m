function [Indicators,Avg_AUC] = FiveFold_SVM(data,label)
%要求数据集每行代表一个样本，label为一行存储了样本标签
Kfold=5;
% M:总样本量； N:一个样本的元素总数
[M,N]=size(data);
% 进行样本的随机分包，indices为具有M个元素的一维向量，第二个参数为label没有问题。
% 每个元素对应一个样本，其中存储了[1,5]的整数，为该样本的分组编号。
indices=crossvalind('Kfold',label,Kfold);

%交叉验证k=5，5个包轮流作为测试集
AUC=zeros(Kfold,1);
for k=1:Kfold
    % 获得test集元素在数据集中对应的单元编号
    test = (indices == k);
    % train集元素的编号为非test元素的编号
    train = ~test;
    
    % 同样地，以下两个函数要求输入的矩阵——行为样本，列为特征。
    model = fitcsvm(data(train,:),label(train,:));
    test_predict = predict(model,data(test,:));
    %计算正确率
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
