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
%����صľ�����Ч�����۲��ֵ
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

%û��������
NoNoiseClass=setdiff(unique(class),-1);
NoNoiseClassNum=size(NoNoiseClass,1);

%����ĳ�ص�TP������FN������FP����
TrueSampleNum=zeros(NoNoiseClassNum,1);
PositveSampleNum=zeros(NoNoiseClassNum,1);
TruePositive=zeros(NoNoiseClassNum,1);
FalseNegative=zeros(NoNoiseClassNum,1);
FalsePositive=zeros(NoNoiseClassNum,1);

%ͳ������
for i=1:NoNoiseClassNum
    TrueSampleNum(i,1)=sum(class==NoNoiseClass(i));
    PositveSampleNum(i,1)=sum(Clustering==NoNoiseClass(i));
    TruePositive(i,1)=sum((class==NoNoiseClass(i))&(Clustering==NoNoiseClass(i)));
    FalseNegative(i,1)=TrueSampleNum(i,1)-TruePositive(i,1);
    FalsePositive(i,1)=PositveSampleNum(i,1)-TruePositive(i,1);
end

%���������Ч��ָ��Accuracy
Indicators.accuracy=sum(TruePositive)/size(class,1);

%���������Ч��ָ��recall
Indicators.recall=TruePositive./(TruePositive+FalseNegative);

%���������Ч��ָ��precision
Indicators.precision=TruePositive./(TruePositive+FalsePositive+eps);

beta=1;
%���������Ч��ָ��F-measure
Indicators.F_measure=((1+beta^2).*Indicators.recall.*Indicators.precision)./...
    (beta^2.*Indicators.recall+Indicators.precision+eps);

%���������Ч��ָ��G-mean
Indicators.G_mean=prod(Indicators.recall);

%���Accuracy
disp(['          Accuracy ',num2str(Indicators.accuracy)]);
%���G-mean
disp(['          G-mean ',num2str(Indicators.G_mean)]);

%������صľ�����Ч����Ϣ
for i=1:NoNoiseClassNum
    
    %����ر��
    disp(['          ',num2str(i),':']);
    
    %���ΪTrue����������
    disp(['              True Samples: ',num2str(TrueSampleNum(i)),...
        ',    Clustered Samples: ',num2str(PositveSampleNum(i))]);
    
    %���ΪTP��FN��FP����������
    disp(['              True Positive: ',num2str(TruePositive(i)),...
        ',    False Negative: ',num2str(FalseNegative(i)),...
        ',    False Positive: ',num2str(FalsePositive(i))]);
    
    %���recall��precision��F-measure��ֵ
    disp(['              recall: ',num2str(Indicators.recall(i)),...
        ',    precision: ',num2str(Indicators.precision(i)),...
        ',    F-measure: ',num2str(Indicators.F_measure(i))]);
    
end

end