function [ chromosome_sub] = sort_sub( weight_sub, chromosomes ,idealp,V,M)
%SORT_SUB Summary of this function goes here
% N=2��������Ⱦɫ�������
[N, ~]= size(chromosomes);
% ��ȡ����Ⱦɫ�����Ӧ�Ⱥ���ֵ
objectives = chromosomes(:,V+1:V+M);
% part2������Ŀǰ����Ӧ�Ⱥ���ֵ��ȫ����С��Ӧ�Ⱥ���ֵ�ľ���ֵ����
part2 = abs(objectives-idealp(ones(N,1),:));
% ����ÿ�е����ֵ����
sub_objectives = max(weight_sub(ones(N,1),:).*part2,[],2);
% ����ÿ��Ⱦɫ����������Ӧ�Ⱥ���ֵ����С���Ǹ�����������Ϊ��Ⱦɫ����������ݡ�
[sorted_obj index] = sort(sub_objectives);
% ��Ⱦɫ���������
chromosome_sub = chromosomes(index(1),:);
end

