function Indicators=EvaluationCriteria_KnownClass(class,Clustering)

%% Evaluation criteria while classes of samples are known
% Input
%       class: given, known classes of samples
%       Clustering: the clustered classes of samples
%       ClusterCenters: The identified cluster centers
% Output
%       Indicators: several evaluation criteria while known classes
%       Indicators.accuracy
%       Indicators.recall
%       Indicators.precision
%       Indicators.F_measure
%       Indicators.G_mean
% Author
%       Junliang Shang (shangjunliang110@163.com)
% Date
%       2017/2/13

%%
%输出簇的聚类有效性评价测度值
disp('###  Evaluation criteria while classes are known  ###');

%%
% order class and clustering which begin from 1, and -1 represents noise.
% ClusterCentersNum=size(ClusterCenters,2);
% 
% class=10000*class;
% class(class<0)=-1;
% for i=1:ClusterCentersNum
%     if class(ClusterCenters(i))<i
%         OldSum=sum((class==class(ClusterCenters(i))) & (Clustering==class(ClusterCenters(i))));
%         NewSum=sum((class==class(ClusterCenters(i))) & (Clustering==i));
%         
%         if NewSum>OldSum
%             class(class==class(ClusterCenters(i)))=i;
%         end
%         
%     else
%         class(class==class(ClusterCenters(i)))=i;
%     end
% end

%没有噪声的
NoNoiseClass=setdiff(unique(class),-1);
NoNoiseClassNum=size(NoNoiseClass,1);

%计算某簇的TP个数，FN个数和FP个数
TrueSampleNum=zeros(NoNoiseClassNum,1);
PositveSampleNum=zeros(NoNoiseClassNum,1);
TruePositive=zeros(NoNoiseClassNum,1);
FalseNegative=zeros(NoNoiseClassNum,1);
FalsePositive=zeros(NoNoiseClassNum,1);

%统计以上
for i=1:NoNoiseClassNum
    TrueSampleNum(i,1)=sum(class==NoNoiseClass(i));
    PositveSampleNum(i,1)=sum(Clustering==NoNoiseClass(i));
    TruePositive(i,1)=sum((class==NoNoiseClass(i))&(Clustering==NoNoiseClass(i)));
    FalseNegative(i,1)=TrueSampleNum(i,1)-TruePositive(i,1);
    FalsePositive(i,1)=PositveSampleNum(i,1)-TruePositive(i,1);
end

%计算聚类有效性指标Accuracy
Indicators.accuracy=sum(TruePositive)/size(class,1);

%计算聚类有效性指标recall
Indicators.recall=TruePositive./(TruePositive+FalseNegative);

%计算聚类有效性指标precision
Indicators.precision=TruePositive./(TruePositive+FalsePositive+eps);

beta=1;
%计算聚类有效性指标F-measure
Indicators.F_measure=((1+beta^2).*Indicators.recall.*Indicators.precision)./...
    (beta^2.*Indicators.recall+Indicators.precision+eps);

%计算聚类有效性指标G-mean
Indicators.G_mean=prod(Indicators.recall);

%输出Accuracy
disp(['          Accuracy ',num2str(Indicators.accuracy)]);
%输出G-mean
disp(['          G-mean ',num2str(Indicators.G_mean)]);

%输出各簇的聚类有效性信息
for i=1:NoNoiseClassNum
    
    %输出簇标号
    disp(['          ',num2str(i),':']);
    
    %输出为True的样本个数
    disp(['              True Samples: ',num2str(TrueSampleNum(i)),...
        ',    Clustered Samples: ',num2str(PositveSampleNum(i))]);
    
    %输出为TP、FN和FP的样本个数
    disp(['              True Positive: ',num2str(TruePositive(i)),...
        ',    False Negative: ',num2str(FalseNegative(i)),...
        ',    False Positive: ',num2str(FalsePositive(i))]);
    
    %输出recall，precision和F-measure的值
    disp(['              recall: ',num2str(Indicators.recall(i)),...
        ',    precision: ',num2str(Indicators.precision(i)),...
        ',    F-measure: ',num2str(Indicators.F_measure(i))]);
    
end

end