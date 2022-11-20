function [over_node,Community] = localexpansion(adj,k)
% 算法的主要框架
% community：cell结构
% 保存重叠的节点
over_node=[];
tic;
NoCore=[];
vertexNum  = size(adj,1);
unassigned = true(1,vertexNum);
Community  = {};
% 原参数 k=3,处理后仍为 k=3；如果原参数k>3,处理后为k，不变
k = max(min(k,vertexNum),3);
% 算法第一步，从网络中随机选择一个点v，通过这个点寻找出该点所在社团的core
% 策略：不断的从当前节点的邻居中找出比当前节点聚类系数大的邻居节点，然后以改邻居点作为新的当前节点；
% 直到不存在比当前节点聚类系数更大的邻居;
% calculate the clustering coefficient of each node
cc = zeros(1,vertexNum);
for i = 1 : vertexNum
    cc(i) = cal_cc(adj,i);%%计算每个节点的聚类系数（度）
end

% now start clutering
% 这个是整个程序的大循环，下面写的是怎样找到core点，找到core后进行聚类操作，直至跳出循环，把已经聚类的节点的标签设置为1，表示已经遍历，
while any(unassigned)
    % 然后从标签为0的节点里再随机找一点进行类似操作，找core聚类
    clc;fprintf('聚类已完成%6.2f%%,耗时%.4f秒... ...\n',sum(~unassigned)/vertexNum*100,toc);
    % find the core node among all the unassigned nodes
    unassignedSet = find(unassigned);
    % 随机从未划分的节点中选择一个作为迭代开始的点
    randIndex     = randi(length(unassignedSet));
    % 记录当前节点的标号
    currentNode   = unassignedSet(randIndex);

    %---Modified by Tian, 7/27/2015---
    %
    %
    
    % 使用当前节点初始化社团的待添加节点集
    addNodes  = currentNode;
    % 找到候选节点，这些节点与当前节点（current node）相连，度大于当前节点，且未被划分
    candidate = find(adj(currentNode,:) & cc > cc(currentNode) & unassigned);
    % 如果没有符合条件的候选节点，则不执行；
    % 否则，......
    while ~isempty(candidate)
        intimateValue = zeros(1,length(candidate));
        % 计算每个候选节点与当前节点的亲密程度
        for i = 1 : length(candidate)
            intimateValue(i) = cal_intimate(adj,candidate(i),addNodes);
        end
        % 选取亲密度最大的节点
        [~,index] = max(intimateValue);
        ci = candidate(index);
        % 将亲密度最大的节点存入社团的待添加节点集中
        addNodes = [addNodes,ci];
        % 将亲密度最大的节点作为新的当前节点，重新选择候选节点
        candidate = find(adj(ci,:) & cc > cc(ci) & unassigned);
    end
    
    % 选择度最大的节点作为核心节点
    [~,max_index] = max(cc(addNodes));
    coreNode =addNodes(max_index);
    
    %---Modified by Tian, 7/27/2015---
    % do local expansion on coreNode
    % 寻找k阶完全图，subCommunity即为我们要找的局部社团（单个）
    subCommunity = find_k_complete(adj,coreNode,k);
    
    % 找到subCommunity中节点的所有邻居
    U = find_neighbors(adj,subCommunity);
    for i=1:length(U)
        a=U(i);
        % intersect函数用法如下：[c, ia, ib] = intersect(A, B);
        % 这个函数是c返回A、B的交集，ia，ib返回的是交集所在数组的指标。
        % 寻找a的邻居中属于subCommunity的节点个数
        edges_a=length(intersect(find(adj(a,:)),subCommunity));
        %%保证强社团,将有60%的边与subCommunity连接的节点加入到局部社团中
        if edges_a/cc(a)>0.6
            subCommunity=[subCommunity a];
        end
    end
    % 评价subCommunity中所有节点，内部连接小于0.5的节点将被加入subCommunityCP
    subCommunityCP = [];
    for i=1:length(subCommunity)
        b=subCommunity(i);
        % C=setdiff(A,B)函数返回在向量A中却不在向量B中的元素，并且C中不包含重复元素，并且从小到大排序。
        % 得到删掉b的subCommunity集，并升序排序
        subCommunity_b=setdiff(subCommunity,b);
        if length(subCommunity_b)==0
            break;
        end
        % 得到b在subCommunity中的邻居
        edges_b=length(intersect(find(adj(b,:)),subCommunity_b));
        % 如果b与subCommunity中节点相连的边仅占其所有边的不到一半，将b节点加入subCommunityCP
        if edges_b/cc(b)<0.5
            subCommunityCP=[subCommunityCP b];
        end
    end
    % 从subCommunity中减去subCommunityCP，并将所得集合中的元素升序排序
    subCommunity=setdiff(subCommunity,subCommunityCP);
    subCommunity=unique([subCommunity coreNode]);
    % subCommunity=intersect(subCommunity,find(unassigned));
    
    % unassigned初始化为全1向量，表示节点没有被划分；
    % 如果subCommunity中包含了已经划分过的点
    if   length(find(unassigned(subCommunity)==0))>0
        % 得到subCommunity中已经划分过的点的标号
        index=find(unassigned(subCommunity)==0);
        % over_node中保存了被多次划分的节点（重叠节点）
        over_node=[over_node subCommunity(index)];

    end
    
    % 将subCommunity中的节点标记为已划分
    unassigned(subCommunity) = false;
    % 将我们找到的一个局部社团加入集合中
    Community = [Community,subCommunity];
    % 生成一个逻辑0的向量，长度与局部社团个数相同
    del = false(1,length(Community));
    
    % Community{end}为Community向量的最后一个元素
    % 除最后一个元素外，元素i与end进行比较，元素i中存在但元素end中不存在的记为1，否则记为0,元素end为新加入的局部社团
    % all函数：检测矩阵中是否全为非零元素，如果全为非零，则返回1
    %% 如果新加入的局部社团与之前加入的某个局部社团完全一致，则删掉之前加入的那个局部社团
    for i = 1 : length(Community) - 1
        del(i) = all(ismember(Community{i},Community{end}));
    end
    Community(del) = [];
    
% end of the whole "while" for discovering all local communities
end

% end of the function localexpansion
end

function cc = cal_cc(adj,node)
% calculate the clustering coefficient of node
%     neighborSet = find(adj(node,:));
%     neighborNum = length(neighborSet);
%     cc = length(find(adj(neighborSet,neighborSet)))/(neighborNum*(neighborNum-1));
%     if isnan(cc)
%         cc = 0;
%     end
%---Modified by Tian, 7/27/2015---
% here just calculate the degree of each node
Binary_node_adj=adj(node,:)~=0;
cc = sum(Binary_node_adj);
%---------------------------------
end

%---Modified by Tian, 7/27/2015---
function value = cal_intimate(adj,A,B)
% calculate the intimate between two sets
common = intersect(find_neighbors(adj,A),find_neighbors(adj,B));
value  = sum(sum(adj(common,common))) ./ length(common) ./ (length(common)-1);
end
%---------------------------------

%---Modified by Tian, 7/27/2015---
function neighbors = find_neighbors(adj,A)
% find all the neighbors of A in adj
% ismember(a,b)依次判断a中元素是不是b中的成员
% 
A = ismember(1:size(adj,1),A);
Binary_node_adj=adj(A,:)~=0;
neighbors = find(any(Binary_node_adj,1) & ~A);
end
%---------------------------------

function neighborSet = find_k_complete(adj,node,k)
% find one k-order complete subgraph in adj which contains node randomly
% 原始k=3，返回所有能与node构成三阶完全图,且除node之外还有至少一个共同邻居的节点，包括node本身

%     neighborSet = [];
%---Modified by Tian, 7/27/2015---
neighborSet = false(1,size(adj,1));
%---------------------------------
% 如果node的邻居个数大于1(这里的node是核心节点)
if length(find(adj(node,:))) > 1
    % node的邻居的所有二阶组合
    allSubMap = nchoosek(find(adj(node,:)),k-1);
    % 将所有二阶组合的顺序打乱
    allSubMap = allSubMap(randperm(size(allSubMap,1)),:);
    % 对于所有二阶组合，
    for i = 1 : size(allSubMap,1)
        % 将核心节点加入组合
        nodes = [node,allSubMap(i, :)];
        % 如果这三个节点构成完全图
        if adj(nodes,nodes) + eye(k)
            % 得到二阶组合的邻居的交集
            common1 = intersect(find_neighbors(adj,nodes(2)),find_neighbors(adj,nodes(3)));
            if length(common1)>1
                %                 neighborSet = allSubMap(i,:);
                %                 break;
                %---Modified by Tian, 7/27/2015---
                % 用true标记该二阶组合的原始标号
                neighborSet(allSubMap(i,:)) = true;
                %---------------------------------
            end
        end
    end
end
%     neighborSet = [neighborSet,node];
%---Modified by Tian, 7/27/2015---
neighborSet = [find(neighborSet),node];
%---------------------------------
end

%---Modified by Tian, 7/27/2015---
% function R = cal_R(adj,A)
% % calculate the R value of A in adj
%     B = any(adj(A,~ismember(1:size(adj,1),A)),2);
%     B = A(B);
%     B = ismember(1:size(adj,1),B);
%     A = ismember(1:size(adj,1),A);
%     R = (sum(sum(adj(B,B)))/2 + sum(sum(adj(B,A & ~B)))) ./ (sum(sum(adj(B,B)))/2 + sum(sum(adj(B,~B))));
% end
%---------------------------------

% function f = cal_f(adj,varargin)
% % calculate the f value of the community
%     a=1;
%     nodeSet = unique(cell2mat(varargin));
%     kin  = sum(sum(adj(nodeSet,nodeSet)));
%     kout = sum(sum(adj(nodeSet,setdiff(1:size(adj,1),nodeSet))));
%     f    = kin/(kin+kout)^a;
% end

% function Q=Calculate_Q(adj,varargin)
% global degree;
% sum=0;
% nodeSet = unique(cell2mat(varargin));
% edges_num = length(find(adj(nodeSet,nodeSet)))/2;
% % edges_num = sum(sum(adj(nodeSet,nodeSet)))/2;
% n=length(nodeSet);
%     for j=1:n
%         a1=nodeSet(j);
%         for k=1:n
%             a2=nodeSet(k);
%             x=degree(a1)*degree(a2)/(2*edges_num);
%             y=adj(a1,a2)-x;
%             sum=sum+y;
%         end;
%     end;
%     Q=sum/(2*edges_num);
% end

