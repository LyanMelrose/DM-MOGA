function AUC = plot_roc( predict, ground_truth )
%  predict       - 分类器对测试集的分类结果
%  ground_truth - 测试集的正确标签,这里只考虑二分类，即0和1
%  auc            - 返回ROC曲线的曲线下的面积

% 计算出ground_truth中正样本的数目pos_num和负样本的数目neg_num
pos_num = sum(ground_truth==1);
neg_num = sum(ground_truth==0);

% 测试集中样本的个数
m=size(ground_truth,1);
% 升序排序
[pre,Index]=sort(predict);
% 将ground_truth继续调整为元素与升序后的predict一一对应。
ground_truth=ground_truth(Index);
% 初始化x和y
x=zeros(m+1,1);
y=zeros(m+1,1);
AUC=0;
% 将初始点设为（1.0, 1.0）
x(1)=1;y(1)=1;

% 依次取出一个样本，计算ROC曲线和AUC值。
for i=2:m
TP=sum(ground_truth(i:m)==1);FP=sum(ground_truth(i:m)==0);
x(i)=FP/neg_num;
y(i)=TP/pos_num;
AUC=AUC+(y(i)+y(i-1))*(x(i-1)-x(i))/2;
end
% 最终点被设为(0,0)
x(m+1)=0;y(m+1)=0;
% 计算AUC值
AUC=AUC+y(m)*x(m)/2;

% 作图
plot(x,y,'lineWidth',1.5);
xlabel('False Positive Rate');
ylabel('True Positive Rate');
end
