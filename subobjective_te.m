function [ obj ] = subobjective_te( weight, objectives, idealpoint )
%SUBOBJECTIVE 
% 多目标分解时切片数为 40，因此 s 为 40。
s = size(weight, 1);
% objsize 保存个体的个数
objsize = size(objectives,1);
% 将 weight 中的 0 转为 0.00001。
weight((weight == 0))=0.00001;
% 个体的个数必须必须与 weight 的长度相同，或等于 1 
if objsize==s 
    part2 = abs(objectives-idealpoint(ones(objsize,1),:));
%     part2 = objectives;
    obj = max(weight.*part2,[],2);
elseif objsize ==1
    % idealpoint 就是 ideap，即全局适应度函数最小值，计算当前个体的适应度值与 idealp 间的绝对值距离
    part2 = abs(objectives-idealpoint);
%  part2 = objectives;
    obj = max(weight.*part2(ones(s, 1),:),[],2);
else
    error('individual size must be same as weight size, or equals 1');
end


end

