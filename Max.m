function com_max=Max(neighbors,a)
% a 为所有节点的所属模块标号
b=[];
% 对于每个邻居
for i=1:length(neighbors)
    % 提取第 i 个 neighbor 其所属模块标号，保存至向量 b 的第 i 个位置。
    b(i)=a(neighbors(i));
end
% e 中保存最多邻居所属的模块
[~,~,e]=mode(b);
% 如果有多个模块，则全部保存在 E 中。然后在其中随机选择一个返回，即 com_max。
E=e{1};
com_max=E(randi(length(E)));


% if length(E)>1
%     com_max=a(k);
% else
%     com_max=E(randi(length(E)));
% end
% end