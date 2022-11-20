clc; clear all;
sumpop=[];
sumclu_assign=[];
summaxnmi=[];
summodnmi=[];
summaxmodnmi=[];
summaxdnmi=[];
% for circle=1:30
for circle=1:10
tic
%% Load network data represented in adjacency matrix.
%% Load network data represented in adjacency matrix.
%%%ziji1 数据，看分层
% load('Datasets\adj_mat_ziji1');
%%%ziji1 数据，看分层125
% load('Datasets\adj_mat_ziji256');
%%%%%%%%%%论文后面的
% load('Datasets\data_test_FB_real');
% load('Datasets\data_test_5clusters');
% Karate Club
% load('Datasets\adjmat_karate_undirected.mat');


% Dolphins
%  load('Datasets\adjmat_dolphins_undirected.mat');

% Football
% load('Datasets\adjmat_football_undirected.mat');

% Polbooks
% load('Datasets\adjmat_polbooks_undirected.mat');

% GNbenchmark
% load('GNBenchmark\GNBenchmark_kout5_1.mat');

% % Extension of GNbenchmark
% network=importdata('GNBenchmark_Ext\network-0.40.dat');
% community=importdata('GNBenchmark_Ext\community-0.40.dat');
% [adj_mat clu_assignment_real] = convert_to_adjmat(community,network);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Karate Club
% load('Datasets\adjmat_karate_undirected.mat');

% Dolphin
% network=importdata('Datasets\dolphin.dat');
%community=importdata('Datasets\dolphin_com.dat');
%[adj_mat clu_assignment_real]= convert_dolphin_to_adjmat(network,community);
%clear network community;

% GNbenchmark
% load('GNBenchmark\GNBenchmark_kout8_0.mat');

% % Extension of GNbenchmark
%  network=importdata('GNBenchmark_Ext\network-0.50.dat');
% community=importdata('GNBenchmark_Ext\community-0.50.dat');
% [adj_mat clu_assignment_real] = convert_extgn_to_adjmat(community,network);

% Karate Club
%  load('Datasets\adjmat_karate_undirected.mat');

% Dolphins
%  load('Datasets\adjmat_dolphins_undirected.mat');

% Football
%  load('Datasets\adjmat_football_undirected.mat');

% Polbooks
% load('Datasets\adjmat_polbooks_undirected.mat');

% GNbenchmark
%load('GNBenchmark\GNBenchmark_kout5_1.mat');

% Extension of GNbenchmark
%network=importdata('GNBenchmark_Ext\network-0.20.dat');
%community=importdata('GNBenchmark_Ext\community-0.20.dat');
%[adj_mat clu_assignment_real] = convert_to_adjmat(community,network);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Karate Club
% load('Datasets\adjmat_karate_undirected.mat');

% Dolphin
% network=importdata('Datasets\dolphin.dat');
%community=importdata('Datasets\dolphin_com.dat');
%[adj_mat clu_assignment_real]= convert_dolphin_to_adjmat(network,community);
%clear network community;

% GNbenchmark
% load('GNBenchmark\GNBenchmark_kout5_0.mat');

% Extension of GNbenchmark
% network=importdata('GNBenchmark_Ext\network-0.50.dat');
% community=importdata('GNBenchmark_Ext\community-0.50.dat');
% [adj_mat clu_assignment_real] = convert_extgn_to_adjmat(community,network);

% adj_mat=load('Datasets\SFI1.txt');
% adj_mat=load('Datasets\netscience_remove1.txt');


% adj_mat = load('E:\我搞定的程序\SI\DPSOCD\email.txt');
% adj_mat = load('E:\我搞定的程序\SI\DPSOCD\PGP_adj.txt');

 X = load('.\数据\footballAdjMatrix');
 adj_mat=X.footballAdjMatrix;% karate为一个结构体，不能直接赋值给AdjMatrix
 clu_assignment_real=[7 0 2 3 7 3 2 8 8 7 3 10 6 2 6 2 7 9 6 1 9 8 8 7 10 0 6 9 11 1 1 6 2 0 6 1 5 0 6 2 3 7 5 6 4  0 11 2 4 11 10 8 3 11 6 1 9 4 11 10 2 6 9 10 2 9 4 11 8 10  9 6 3 11 3 4 9 8 8 1 5 3 5 11 3 6 4 9 11 0 5 4 4 7 1 9 9 10 3 6 2 1 3 0 7 0 2 3 8 0 4 8 4 9 11]+1;

%adj_mat=Adjreverse(6,0);
%adj_mat=[14 4 1 1;4 4 0 1;1 0 16 3;1 1 3 2];
%clu_assignment_real(1:32)=ones(1,32);
%clu_assignment_real(33:64)=ones(1,32)+1;
%clu_assignment_real(65:96)=ones(1,32)+2;
%clu_assignment_real(97:128)=ones(1,32)+3;

starttime = clock;

[V,V] = size(adj_mat); 
edgeslist = edges_list(adj_mat,V);
M = 2;
popsize = 100;
niche =40;
gmax=100;
% gmax=50;
pc=1;
% pc=1;
% pm=0.05;
pm=0.1;
pg=0.9;    %%%%%%%%%%看是从整个领域还是整个种群当中选解
%%%%%%%%%%%%%%%%%%%%%%
% kk=2;
kk=2;
%%%%%%%%%%%%%%%%%%%%
ll=0;
for i=1:V
    for j=i:V
        ll=adj_mat(i,j)+ll;
    end
end
global idealp weights neighbors;
idealp = -Inf*ones(1,M);

[weights neighbors] = init_weight(popsize, niche);
popsize = size(weights,1);
chromosomes = initialize_variables(popsize,edgeslist,adj_mat,M,V,ll);  % 初始化种群

for icounter = 1:gmax
for i=1:popsize
    i_neighbor_index = neighbors(i,:);
    i_neighbor_weight = weights(i_neighbor_index,:);
    i_neighbor_chromosome = chromosomes(i_neighbor_index,:);
    
    old_objective = i_neighbor_chromosome(:,V+1:V+M);
    new_chromosomes = genetic_op(i_neighbor_chromosome,chromosomes,adj_mat,edgeslist,pc,pm,V,M,ll);
    %update idealp
    idealp = min(min(new_chromosomes(:,V+1:V+M)),idealp);
    %update the neighbours
    %选择对于当前子问题最好的解
    chromosome_sub = sort_sub(weights(i,:),new_chromosomes,idealp,V,M);
    %-------------------------
    new_subobj = subobjective_te(i_neighbor_weight, chromosome_sub(V+1:V+M), idealp);
    old_subobj = subobjective_te(i_neighbor_weight, old_objective, idealp);
    jj=0;
    iposition = find(new_subobj<old_subobj);
    jj=length(iposition);
    if jj>kk;                 %%挑取2个进行替代
        kk1=floor(jj*rand(kk,1))+1;%%rand（2，1）产生随机的2行一列随机数
        iposition1=iposition(kk1);
    else
        iposition1=iposition;    
    end
    
    chromosomes(i_neighbor_index(iposition1)',:) = chromosome_sub(ones(length(iposition1),1),:);   
end
% fprintf('iteration %u finished, time used: %u\n', icounter, toc);
for i = 1:popsize
    clu_assign(i).p = decode(chromosomes(i,1:V));
    modnmi(1,i) = modularity(adj_mat,clu_assign(i).p);
    modnmi(2,i) = normalized_mutual_information(clu_assignment_real,clu_assign(i).p);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    modnmi(3,i)=modularity_density(clu_assign(i).p,adj_mat,i);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
[max_modularity,max_index1] = max(modnmi(1,:));
[max_nmi,max_index2] = max(modnmi(2,:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[max_modularity_density,max_index3] = max(modnmi(3,:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fprintf('%d-the max_modularity: %u , the correspond nmi: %u\n', circle, max_modularity, modnmi(2,max_index1));
% %%%%%%%%%%%%%%%%%%%%%%%%
% fprintf('%d-the max_modularity_density: %u\n',circle, max_modularity_density);

% fprintf('%d-the max_modularity_density: %u , the correspond nmi: %u\n',circle, max_modularity_density, modnmi(2,max_index3));
%%%%%%%%%%%%%%%%%%%%%%%%%

% fprintf('%d-the max_nmi: %u , the correspond modularity: %u\n',circle, max_nmi, modnmi(1,max_index2));
end

%-------------------------------------
%删除chromosomes中相同的解
objectives = chromosomes(:,V+1:V+M);
Ns = popsize;
% i=1;
% while i<Ns
%    deltf = objectives - ones(Ns,1)*objectives(i,:);
%    deltf(i,:) = inf;
%    aa = find(sum(abs(deltf),2)==0);
%    objectives(aa,:) = []; chromosomes(aa,:) = [];
%    [Ns,temp] = size(objectives); i=i+1;
% end
% -------------------------------------

%Find the solution having the maximum modularity
for i = 1:Ns
    clu_assign(i).p = decode(chromosomes(i,1:V));
    modnmi(1,i) = modularity(adj_mat,clu_assign(i).p);
    modnmi(2,i) = normalized_mutual_information(clu_assignment_real,clu_assign(i).p);
    modnmi(3,i)=modularity_density(clu_assign(i).p,adj_mat,i);
end
[max_modularity,max_index1] = max(modnmi(1,:));
[max_nmi,max_index2] = max(modnmi(2,:));
[max_modularity_density,max_index3] = max(modnmi(3,:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aaa=[];
bbb=[];
aaa1=[];
bbb1=[];
aaa2=[];
bbb2=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[aaa1,bbb1]=sort(modnmi(1,:));
[aaa,bbb]=sort(modnmi(2,:));
[aaa2,bbb2]=sort(modnmi(3,:));%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%bbb是值的顺序nmi，看是否有层次结构
final_clu_assignment1 = clu_assign(max_index1).p;
final_clu_assignment2 = clu_assign(max_index2).p;
final_clu_num1 = max(final_clu_assignment1);
final_clu_num2 = max(final_clu_assignment2);

fprintf('the final solution having the maximum modularity:\n');
for i = 1:final_clu_num1
    fprintf('Cluster %d:',i);
    display(find(final_clu_assignment1 == i));
end
fprintf('objective1: %f; objective2: %f;\nmodularity:%f; nmi:%f;\n',...
chromosomes(max_index1,V+1),chromosomes(max_index1,V+2),...
modnmi(1,max_index1),modnmi(2,max_index1));

fprintf('the final solution having the maximum nmi:\n');
for i = 1:final_clu_num2
    fprintf('Cluster %d:',i);
    display(find(final_clu_assignment2 == i));
end
fprintf('objective1: %f; objective2: %f;\nmodularity:%f; nmi:%f;\n',...
chromosomes(max_index2,V+1),chromosomes(max_index2,V+2),...
modnmi(1,max_index2),modnmi(2,max_index2));

sumpop=[sumpop;chromosomes];
sumclu_assign=[sumclu_assign;clu_assign];
summodnmi=[summodnmi;modnmi];
summaxmodnmi=[summaxmodnmi; max_modularity,modnmi(2,max_index1)]; 
summaxdnmi=[summaxdnmi; max_modularity_density,modnmi(2,max_index3)]; 
summaxnmi=[summaxnmi,max_nmi]

mod_circle(circle)=max_modularity
nmi_circle(circle)=max_nmi
end
savedata1('result60.txt',[mod_circle' nmi_circle']); 
pp=chromosomes(:,V+1:V+M);
plot(pp(:,1)',pp(:,2)','r.');
xlabel('Negative Ratio Association');
ylabel('Ratio Cut');
set(get(gca,'xlabel'),'fontsize',12);
set(get(gca,'ylabel'),'fontsize',12);
modstr = num2str(modnmi(2,:)');
%  text(pp(:,1)',pp(:,2)',modstr);
disp(sprintf('total time used %u', etime(clock, starttime)));
