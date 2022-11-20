clear;
clc;
%%
%
%加载数据
data=load('C:\Desktop\DeepLearningFloder\IrisSet.txt');  
label = [zeros(50,1);ones(50,1)];
%%
train_data = [data(1:25,1:2);data(51:75,1:2)];
train_label = [label(1:25,1);label(51:75,1)];
test_data = [data(26:50,1:2);data(76:100,1:2)];
test_label = [label(26:50,1);label(76:100,1)];
 
%% 建立模型并预测  采用matlab自带的
SVMSVMStruct=svmtrain(train_data,train_label,'showplot',true)    %注意matlab自带的svm可视化只能2-D显示
Group=svmclassify(SVMStruct,test_data,'showplot',true);
