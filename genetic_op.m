function [ f ] = genetic_op( parent_chromosome,chromosomes,edgeslist,pc,pm,V,M,degree1,matrix,Node )
% function [ f ] = genetic_op( parent_chromosome,chromosomes,edgeslist,pc,pm,V,M,ll,degree1,matrix,Node )
%GENETIC_OP Summary of this function goes here

[N,~] = size(parent_chromosome);
[N1,~]=size(chromosomes);
r0=0;


% crossover（交叉算子）
p = 1;
for i = 1 : 1
    % Select the first parent
    parent_1 = ceil(N*rand(1));
    r0=rand(1);
    if r0<=0.9
        % Select the second parent
        parent_2 = ceil(N*rand(1));
        % Make sure both the parents are not the same.
        while parent_2 == parent_1
            parent_2 = ceil(N*rand(1));
        end
        parent_1 = parent_chromosome(parent_1,:);
        parent_2 = parent_chromosome(parent_2,:);
    else
        parent_2=N+floor((N1-N)*rand(1));
        parent_1 = parent_chromosome(parent_1,:);
        parent_2 = chromosomes(parent_2,:);
    end
    % Get the chromosome information for each randomnly selected parents
    % parent_1 = parent_chromosome(parent_1,:);
    % parent_2 = parent_chromosome(parent_2,:);
    if rand<=pc
        % Perform corssover for each decision variable in the chromosome.
        child_1 = parent_1(1:V);
        child_2 = parent_2(1:V);
        position1 = ceil(V*rand(1));
        position2 = ceil(V*rand(1));
        while position1==position2
            position2 = ceil(V*rand(1));
        end
        %         while position2 == position1
        %             position2 = ceil(V*rand(1));
        %         end
        minp = min(position1,position2);
        maxp = max(position1,position2);
        r1=0;
        r1=floor(3*rand(1)+1);   %%floor（A）取A最小整数
        %%分三段随机交叉
        if r1==1
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            temp = child_1(1:minp);
            child_1(1:minp) = child_2(1:minp);
            child_2(1:minp) = temp;
            % Evaluate the objective function for the offsprings and as before
            % concatenate the offspring chromosome with objective value.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/7/30%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            child_1(:,V + 1: M + V) = evaluate_objective(child_1,degree1,matrix,Node);
            child_2(:,V + 1: M + V) = evaluate_objective(child_2,degree1,matrix,Node);
            %         child_1(:,V + 1: M + V) = evaluate_objective(child_1,ll,degree1,matrix,Node);
            %         child_2(:,V + 1: M + V) = evaluate_objective(child_2,ll,degree1,matrix,Node);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif r1==2
            temp = child_1(minp:maxp);
            child_1(minp:maxp) = child_2(minp:maxp);
            child_2(minp:maxp) = temp;
            % Evaluate the objective function for the offsprings and as before
            % concatenate the offspring chromosome with objective value.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/7/30%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            child_1(:,V + 1: M + V) = evaluate_objective(child_1,degree1,matrix,Node);
            child_2(:,V + 1: M + V) = evaluate_objective(child_2,degree1,matrix,Node);
            %             child_1(:,V + 1: M + V) = evaluate_objective(child_1,ll,degree1,matrix,Node);
            %             child_2(:,V + 1: M + V) = evaluate_objective(child_2,ll,degree1,matrix,Node);
        elseif r1==3
            temp = child_1(maxp:V);
            child_1(maxp:V) = child_2(maxp:V);
            child_2(maxp:V) = temp;
            % Evaluate the objective function for the offsprings and as before
            % concatenate the offspring chromosome with objective value.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/7/30%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            child_1(:,V + 1: M + V) = evaluate_objective(child_1,degree1,matrix,Node);
            child_2(:,V + 1: M + V) = evaluate_objective(child_2,degree1,matrix,Node);
            %             child_1(:,V + 1: M + V) = evaluate_objective(child_1,ll,degree1,matrix,Node);
            %             child_2(:,V + 1: M + V) = evaluate_objective(child_2,ll,degree1,matrix,Node);
        end
    else
        child_1 = parent_1(1:V+M);
        child_2 = parent_2(1:V+M);
    end
    f(p,:) = child_1;
    f(p+1,:) = child_2;
    p = p + 2;
end

% mutation（变异算子）
[NN,~] = size(f);
for i = 1 : NN
    for j = 1:V
        if rand<=0.1
            if edgeslist(j).n >= 2  %%邻接点有2个以上
                temp = f(i,j);
                while 1
                    f(i,j) =edgeslist(j).e(ceil(rand*edgeslist(j).n));
                    if temp ~= f(i,j)
                        break;
                    end
                end
            end
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/7/30%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    f(i,V + 1: M + V) = evaluate_objective(f(i,1:V),degree1,matrix,Node);
    %     f(i,V + 1: M + V) = evaluate_objective(f(i,1:V),ll,degree1,matrix,Node);
end
end

