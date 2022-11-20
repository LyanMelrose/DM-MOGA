function [ chromosome_sub] = sort_sub( weight_sub, chromosomes ,idealp,V,M)
%SORT_SUB Summary of this function goes here
% N=2，有两条染色体待排序。
[N, ~]= size(chromosomes);
% 提取两条染色体的适应度函数值
objectives = chromosomes(:,V+1:V+M);
% part2，计算目前的适应度函数值到全局最小适应度函数值的绝对值距离
part2 = abs(objectives-idealp(ones(N,1),:));
% 搜索每行的最大值返回
sub_objectives = max(weight_sub(ones(N,1),:).*part2,[],2);
% 根据每个染色体中两个适应度函数值中最小的那个进行排序，作为对染色体排序的依据。
[sorted_obj index] = sort(sub_objectives);
% 对染色体进行排序。
chromosome_sub = chromosomes(index(1),:);
end

