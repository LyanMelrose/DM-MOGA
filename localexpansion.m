function [over_node,Community] = localexpansion(adj,k)
% �㷨����Ҫ���
% community��cell�ṹ
% �����ص��Ľڵ�
over_node=[];
tic;
NoCore=[];
vertexNum  = size(adj,1);
unassigned = true(1,vertexNum);
Community  = {};
% ԭ���� k=3,�������Ϊ k=3�����ԭ����k>3,�����Ϊk������
k = max(min(k,vertexNum),3);
% �㷨��һ���������������ѡ��һ����v��ͨ�������Ѱ�ҳ��õ��������ŵ�core
% ���ԣ����ϵĴӵ�ǰ�ڵ���ھ����ҳ��ȵ�ǰ�ڵ����ϵ������ھӽڵ㣬Ȼ���Ը��ھӵ���Ϊ�µĵ�ǰ�ڵ㣻
% ֱ�������ڱȵ�ǰ�ڵ����ϵ��������ھ�;
% calculate the clustering coefficient of each node
cc = zeros(1,vertexNum);
for i = 1 : vertexNum
    cc(i) = cal_cc(adj,i);%%����ÿ���ڵ�ľ���ϵ�����ȣ�
end

% now start clutering
% �������������Ĵ�ѭ��������д���������ҵ�core�㣬�ҵ�core����о��������ֱ������ѭ�������Ѿ�����Ľڵ�ı�ǩ����Ϊ1����ʾ�Ѿ�������
while any(unassigned)
    % Ȼ��ӱ�ǩΪ0�Ľڵ����������һ��������Ʋ�������core����
    clc;fprintf('���������%6.2f%%,��ʱ%.4f��... ...\n',sum(~unassigned)/vertexNum*100,toc);
    % find the core node among all the unassigned nodes
    unassignedSet = find(unassigned);
    % �����δ���ֵĽڵ���ѡ��һ����Ϊ������ʼ�ĵ�
    randIndex     = randi(length(unassignedSet));
    % ��¼��ǰ�ڵ�ı��
    currentNode   = unassignedSet(randIndex);

    %---Modified by Tian, 7/27/2015---
    %
    %
    
    % ʹ�õ�ǰ�ڵ��ʼ�����ŵĴ���ӽڵ㼯
    addNodes  = currentNode;
    % �ҵ���ѡ�ڵ㣬��Щ�ڵ��뵱ǰ�ڵ㣨current node���������ȴ��ڵ�ǰ�ڵ㣬��δ������
    candidate = find(adj(currentNode,:) & cc > cc(currentNode) & unassigned);
    % ���û�з��������ĺ�ѡ�ڵ㣬��ִ�У�
    % ����......
    while ~isempty(candidate)
        intimateValue = zeros(1,length(candidate));
        % ����ÿ����ѡ�ڵ��뵱ǰ�ڵ�����̶ܳ�
        for i = 1 : length(candidate)
            intimateValue(i) = cal_intimate(adj,candidate(i),addNodes);
        end
        % ѡȡ���ܶ����Ľڵ�
        [~,index] = max(intimateValue);
        ci = candidate(index);
        % �����ܶ����Ľڵ�������ŵĴ���ӽڵ㼯��
        addNodes = [addNodes,ci];
        % �����ܶ����Ľڵ���Ϊ�µĵ�ǰ�ڵ㣬����ѡ���ѡ�ڵ�
        candidate = find(adj(ci,:) & cc > cc(ci) & unassigned);
    end
    
    % ѡ������Ľڵ���Ϊ���Ľڵ�
    [~,max_index] = max(cc(addNodes));
    coreNode =addNodes(max_index);
    
    %---Modified by Tian, 7/27/2015---
    % do local expansion on coreNode
    % Ѱ��k����ȫͼ��subCommunity��Ϊ����Ҫ�ҵľֲ����ţ�������
    subCommunity = find_k_complete(adj,coreNode,k);
    
    % �ҵ�subCommunity�нڵ�������ھ�
    U = find_neighbors(adj,subCommunity);
    for i=1:length(U)
        a=U(i);
        % intersect�����÷����£�[c, ia, ib] = intersect(A, B);
        % ���������c����A��B�Ľ�����ia��ib���ص��ǽ������������ָ�ꡣ
        % Ѱ��a���ھ�������subCommunity�Ľڵ����
        edges_a=length(intersect(find(adj(a,:)),subCommunity));
        %%��֤ǿ����,����60%�ı���subCommunity���ӵĽڵ���뵽�ֲ�������
        if edges_a/cc(a)>0.6
            subCommunity=[subCommunity a];
        end
    end
    % ����subCommunity�����нڵ㣬�ڲ�����С��0.5�Ľڵ㽫������subCommunityCP
    subCommunityCP = [];
    for i=1:length(subCommunity)
        b=subCommunity(i);
        % C=setdiff(A,B)��������������A��ȴ��������B�е�Ԫ�أ�����C�в������ظ�Ԫ�أ����Ҵ�С��������
        % �õ�ɾ��b��subCommunity��������������
        subCommunity_b=setdiff(subCommunity,b);
        if length(subCommunity_b)==0
            break;
        end
        % �õ�b��subCommunity�е��ھ�
        edges_b=length(intersect(find(adj(b,:)),subCommunity_b));
        % ���b��subCommunity�нڵ������ı߽�ռ�����бߵĲ���һ�룬��b�ڵ����subCommunityCP
        if edges_b/cc(b)<0.5
            subCommunityCP=[subCommunityCP b];
        end
    end
    % ��subCommunity�м�ȥsubCommunityCP���������ü����е�Ԫ����������
    subCommunity=setdiff(subCommunity,subCommunityCP);
    subCommunity=unique([subCommunity coreNode]);
    % subCommunity=intersect(subCommunity,find(unassigned));
    
    % unassigned��ʼ��Ϊȫ1��������ʾ�ڵ�û�б����֣�
    % ���subCommunity�а������Ѿ����ֹ��ĵ�
    if   length(find(unassigned(subCommunity)==0))>0
        % �õ�subCommunity���Ѿ����ֹ��ĵ�ı��
        index=find(unassigned(subCommunity)==0);
        % over_node�б����˱���λ��ֵĽڵ㣨�ص��ڵ㣩
        over_node=[over_node subCommunity(index)];

    end
    
    % ��subCommunity�еĽڵ���Ϊ�ѻ���
    unassigned(subCommunity) = false;
    % �������ҵ���һ���ֲ����ż��뼯����
    Community = [Community,subCommunity];
    % ����һ���߼�0��������������ֲ����Ÿ�����ͬ
    del = false(1,length(Community));
    
    % Community{end}ΪCommunity���������һ��Ԫ��
    % �����һ��Ԫ���⣬Ԫ��i��end���бȽϣ�Ԫ��i�д��ڵ�Ԫ��end�в����ڵļ�Ϊ1�������Ϊ0,Ԫ��endΪ�¼���ľֲ�����
    % all���������������Ƿ�ȫΪ����Ԫ�أ����ȫΪ���㣬�򷵻�1
    %% ����¼���ľֲ�������֮ǰ�����ĳ���ֲ�������ȫһ�£���ɾ��֮ǰ������Ǹ��ֲ�����
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
% ismember(a,b)�����ж�a��Ԫ���ǲ���b�еĳ�Ա
% 
A = ismember(1:size(adj,1),A);
Binary_node_adj=adj(A,:)~=0;
neighbors = find(any(Binary_node_adj,1) & ~A);
end
%---------------------------------

function neighborSet = find_k_complete(adj,node,k)
% find one k-order complete subgraph in adj which contains node randomly
% ԭʼk=3��������������node����������ȫͼ,�ҳ�node֮�⻹������һ����ͬ�ھӵĽڵ㣬����node����

%     neighborSet = [];
%---Modified by Tian, 7/27/2015---
neighborSet = false(1,size(adj,1));
%---------------------------------
% ���node���ھӸ�������1(�����node�Ǻ��Ľڵ�)
if length(find(adj(node,:))) > 1
    % node���ھӵ����ж������
    allSubMap = nchoosek(find(adj(node,:)),k-1);
    % �����ж�����ϵ�˳�����
    allSubMap = allSubMap(randperm(size(allSubMap,1)),:);
    % �������ж�����ϣ�
    for i = 1 : size(allSubMap,1)
        % �����Ľڵ�������
        nodes = [node,allSubMap(i, :)];
        % ����������ڵ㹹����ȫͼ
        if adj(nodes,nodes) + eye(k)
            % �õ�������ϵ��ھӵĽ���
            common1 = intersect(find_neighbors(adj,nodes(2)),find_neighbors(adj,nodes(3)));
            if length(common1)>1
                %                 neighborSet = allSubMap(i,:);
                %                 break;
                %---Modified by Tian, 7/27/2015---
                % ��true��Ǹö�����ϵ�ԭʼ���
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

