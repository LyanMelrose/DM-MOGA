function [ obj ] = subobjective_te( weight, objectives, idealpoint )
%SUBOBJECTIVE 
% ��Ŀ��ֽ�ʱ��Ƭ��Ϊ 40����� s Ϊ 40��
s = size(weight, 1);
% objsize �������ĸ���
objsize = size(objectives,1);
% �� weight �е� 0 תΪ 0.00001��
weight((weight == 0))=0.00001;
% ����ĸ������������ weight �ĳ�����ͬ������� 1 
if objsize==s 
    part2 = abs(objectives-idealpoint(ones(objsize,1),:));
%     part2 = objectives;
    obj = max(weight.*part2,[],2);
elseif objsize ==1
    % idealpoint ���� ideap����ȫ����Ӧ�Ⱥ�����Сֵ�����㵱ǰ�������Ӧ��ֵ�� idealp ��ľ���ֵ����
    part2 = abs(objectives-idealpoint);
%  part2 = objectives;
    obj = max(weight.*part2(ones(s, 1),:),[],2);
else
    error('individual size must be same as weight size, or equals 1');
end


end

