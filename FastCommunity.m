function [Community,MaxDeltaQ] = FastCommunity(WeightNet)
WeightNet = load('RealWorld\dolphin.txt');
Datalabel=load('RealWorld\real_label_dolphin.txt');
% 	if nargin < 1
% 		WeightNet = [1 1/2 1/2 0 0 0;
% 					1/2 3/4 1/4 1/4 1/4 0;
% 					1/2 1/4 3/4 1/4 1/4 0;
% 					0 1/4 1/4 3/4 1/4 1/2;
% 					0 1/4 1/4 1/4 3/4 1/2;
% 					0 0 0 1/2 1/2 1;];
% 	end
	[~,N] = size(WeightNet);
	DiagnolWeight = diag(WeightNet);
	NodeTravelled = zeros(N,1);
	SumWeight = (sum(sum(WeightNet)) + sum(WeightNet(logical(eye(N))))) / 2; 
	Community = {};
    MaxDeltaQ = 0;
	while ~all(NodeTravelled)
		NoSelectNode = find(NodeTravelled == 0);
        numel(NoSelectNode);
		StartNode = NoSelectNode(randi(length(NoSelectNode)));
		NodeTravelled(StartNode) = 1;
		CurrentCommunity = [StartNode];
		while 1
			NeighborNode = [];
			for i = 1 : length(CurrentCommunity)
				NeighborNode = [ NeighborNode find(WeightNet(CurrentCommunity(i),:)~=0)]; 
			end
			ComNodeall = unique(NeighborNode);   %寻找当前Community的邻居节�?包含自循�?
			NeighborNode = setdiff(ComNodeall,CurrentCommunity);   %去掉邻居点中自循�?			
            CommunityNet = WeightNet(CurrentCommunity,CurrentCommunity);   %构�?Community网络
			SigmaInsideWeight = (sum(sum(CommunityNet)) + sum(CommunityNet(logical(eye(length(CurrentCommunity)))))) / 2;
			KNeighborInside = zeros(length(NeighborNode),1);
			KNeighborTotal = zeros(length(NeighborNode),1);
			for i = 1 : length(NeighborNode)
				KNeighborTotal(i) = sum(WeightNet(NeighborNode(i),:));
				KNeighborInside(i) = sum(WeightNet(NeighborNode(i),[CurrentCommunity]));
			end
			SigmaOutsideWeight = sum(KNeighborInside); 
%			CurrentQ = sum(sum(CommunityNet)) / 2 + sum(CommunityNet(logical(eye(length(CommunityNet))))) / 2; 
%			CurrentQ = CurrentQ / (2*N);
            DeltaQ = zeros(length(NeighborNode),1);
			for i = 1 : length(NeighborNode)
 				DeltaQ(i) = (((SigmaInsideWeight + KNeighborInside(i)) / (2*SumWeight)) - ((SigmaInsideWeight + SigmaOutsideWeight + KNeighborTotal(i)) / (2*SumWeight))^2) - ((SigmaInsideWeight / (2*SumWeight)) - ((SigmaInsideWeight + SigmaOutsideWeight) / (2*SumWeight))^2 - (KNeighborTotal(i) / (2*SumWeight))^2);
%				NeighborCommunity = [];
%				NeighborCommunity = [CurrentCommunity NeighborNode(i)];
%				NeighborCommunityNet = WeightNet(NeighborCommunity , NeighborCommunity) - 4/(2*N);
%				NeighborQ = (sum(sum(NeighborCommunityNet)) + sum(NeighborCommunityNet(logical(eye(length(NeighborCommunity))))))/2;
%				NeighborQ = NeighborQ / (2*N);
%                DeltaQ(i) = NeighborQ - CurrentQ;
            end
            
			[CurrentMaxDeltaQ,MaxNode] = max(DeltaQ);
			if CurrentMaxDeltaQ > MaxDeltaQ
                MaxDeltaQ = CurrentMaxDeltaQ;
            end
            if CurrentMaxDeltaQ < 0 
				break;
			end
			CurrentCommunity = [CurrentCommunity,NeighborNode(MaxNode)];
			NodeTravelled(NeighborNode(MaxNode)) = 1;
            if isempty(NeighborNode)
                break;
            end
		end
		Flagbit = 0;
		for i = 1 : length(Community)
			if isequal(intersect(Community{i},CurrentCommunity),Community{i})
				Community(i) = {CurrentCommunity};
				Flagbit = 1;
			end
		end
		if Flagbit == 0
			Community = [Community ; {CurrentCommunity}];
		end
	end
end
