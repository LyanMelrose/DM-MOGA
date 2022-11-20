clear;
clc;
%%
%
%��������
data=load('C:\Desktop\DeepLearningFloder\IrisSet.txt');  
label = [zeros(50,1);ones(50,1)];
%%
train_data = [data(1:25,1:2);data(51:75,1:2)];
train_label = [label(1:25,1);label(51:75,1)];
test_data = [data(26:50,1:2);data(76:100,1:2)];
test_label = [label(26:50,1);label(76:100,1)];
 
%% ����ģ�Ͳ�Ԥ��  ����matlab�Դ���
SVMSVMStruct=svmtrain(train_data,train_label,'showplot',true)    %ע��matlab�Դ���svm���ӻ�ֻ��2-D��ʾ
Group=svmclassify(SVMStruct,test_data,'showplot',true);
