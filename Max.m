function com_max=Max(neighbors,a)
% a Ϊ���нڵ������ģ����
b=[];
% ����ÿ���ھ�
for i=1:length(neighbors)
    % ��ȡ�� i �� neighbor ������ģ���ţ����������� b �ĵ� i ��λ�á�
    b(i)=a(neighbors(i));
end
% e �б�������ھ�������ģ��
[~,~,e]=mode(b);
% ����ж��ģ�飬��ȫ�������� E �С�Ȼ�����������ѡ��һ�����أ��� com_max��
E=e{1};
com_max=E(randi(length(E)));


% if length(E)>1
%     com_max=a(k);
% else
%     com_max=E(randi(length(E)));
% end
% end