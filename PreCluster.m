function [ Community1,index ] = PreCluster( M )
%PreCluster Function Declaration
%   Input:
%       AdjMatrix : The adjacency matrix, 0-1 stands for Non-joined or Joined, respectively. And No self-circulation.
%       Point_Seq : The input sequence of Points
%       Limitnum : The Biggest Clique point Number
%   Output:
%       Community : A partition of the Community. Type: cell

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Pre Processing
%        M= [0 1 1 1 1 1 0;
%            1 0 0 1 0 0 0;
%            1 0 0 1 1 0 0;
%            1 1 1 0 1 0 1;
%            1 0 1 1 0 1 1;
%            1 0 0 0 1 0 1;
%            0 0 0 1 1 1 0;];
       
        AdjMatrix = M;
          
           LimitNum = 4;
           Point_Seq = 1:size(AdjMatrix, 1);
        
 
%         fprintf('Error: Wrong Input!\n');
%         fprintf('Usage: [ Community ] = PreCluster( AdjMatrix, Point_Seq, LimitNum )\n');
%         fprintf('    or [ Community ] = PreCluster(  )\n');
%         return;
  
%    if length(AdjMatrix) ~= length(Point_Seq)
%        fprintf('Error: The length of AdjMatrix and the Length of Point_Seq Must be equal!\n');
%        return;
%    end
    Point_Num = size(Point_Seq, 2);
    if AdjMatrix(1,1) == 1
        AdjMatrix(logical(eye(Point_Num))) = 0;
    end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%start
index=0;
    Degree = sum(AdjMatrix,1);    
    i = 1;
    while i <= Point_Num
        Point_Current_No = Point_Seq(i);
        Point_Current_Neighbours = find(AdjMatrix(Point_Current_No,:));      %Finde Current Point Neighbours Points
        Degree_Current_Neighbours = Degree(Point_Current_Neighbours);
        [~ , Ponit_MaxDeg_No] = max(Degree_Current_Neighbours);              %Find the Biggest degree From neighbours Point
        Point_Next_No = Point_Current_Neighbours(Ponit_MaxDeg_No);
        Community_Current=[Point_Current_No,Point_Next_No];                  %Combine Current Community
        while any( all( AdjMatrix( Community_Current' , : ) , 1 ) ) && length(Community_Current) < LimitNum          %At least For 3 Points Community If length of Community>=3 ,then Save it
            Point_Shared_Neighbours = find( all( AdjMatrix( Community_Current' , : ) , 1 ) );%Continue to extend the Community by Finding shared Neighbours
            Degree_Shared_Neighbours = Degree(Point_Shared_Neighbours);
            [~ , Ponit_MaxDeg_No] = max(Degree_Shared_Neighbours);
            Point_Next_No = Point_Shared_Neighbours(Ponit_MaxDeg_No);
            Community_Current = [Community_Current,Point_Next_No];
        end
        if length(Community_Current) >3
            Community{i} = Community_Current;
            index=1;
        else
%            Community{i} = [Point_Current_No];
        end
        i = i+1;
    end
    if index==0
        Community1=[];
    else
        for j=1:length(Community)
            if length(Community{j})==4
               Community1=Community{j};
            end
        end
    end
    
end

