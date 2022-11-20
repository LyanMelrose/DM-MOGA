clear 
clc
tic
global numObjectives popsize numVar niche max_gen mutate_posibility  edgeslist V  %多目标全局参数
global strNetwork AdjacentMatrix Datalabel  idealp weights neighbors  ll M  adj_mat  pc pg pm             %网络聚类参数
global Thirty_Run_maxQ Thirty_Run_maxNMI      %实验参数
%% 多目标参数赋予初值

max_gen         = 100;       % maximal number of generations
niche           = 40;        % neighborhood size
popsize			= 100;		 % population size
runtimes        = 10;		 % run times
M = 2;
pg=1;
pm=0.9;
pc=1;

%% 网络参数 numVar;邻接矩阵；真实划分；每个点的信息node

strNetwork='karate';
AdjacentMatrix = load('RealWorld\karate.txt');
Datalabel=load('RealWorld\real_label_karate.txt');
adj_mat=AdjacentMatrix;
[V,V] = size(adj_mat); 
edgeslist = edges_list(adj_mat,V);

numVar=size(AdjacentMatrix,1);
 ll=0;
for i=1:V
    for j=i:V
        ll=adj_mat(i,j)+ll;
    end
end

%% 实验部分
root = sprintf('results1/%s/ParetoFront',strNetwork);  %%建个文件保存实验数据
if ~isdir(root) %判断路径是否存在
    mkdir(root); 
end
root = sprintf('results1/%s/metrics',strNetwork);
if ~isdir(root) %判断路径是否存在
    mkdir(root);
end

Thirty_Run_maxQ=[];
Thirty_Run_maxNMI=[];

for mg=1:runtimes
    fprintf('第%s次运行\n',num2str(mg));
    TMOEAD(mg);
    fprintf('\n');
end 

%[~,index1]=max(Thirty_Run_maxQ(:,1));%找到最大的Q索引
%[~,index2]=max(Thirty_Run_maxNMI(:,2));%找到最大的NMI索引

fprintf('Thirtyrun_results:\n');
fprintf('maxQ =    %g %g\n',Thirty_Run_maxQ(index1,1),Thirty_Run_maxQ(index1,2));    %输出最大的Q时的Q和NMI
fprintf('maxNMI =  %g %g\n',Thirty_Run_maxNMI(index2,1),Thirty_Run_maxNMI(index2,2));%输出最大的NMI时的NMI和Q

path = sprintf('results1/%s/metrics/MODPSO1_%s_Thirty_Run_maxQ.txt',strNetwork,strNetwork);%保存实验数据取平均值
savedata1(path,[Thirty_Run_maxQ;0 0;mean(Thirty_Run_maxQ(:,1)) mean(Thirty_Run_maxQ(:,2))]); 
path = sprintf('results1/%s/metrics/MODPSO1_%s_Thirty_Run_maxNMI.txt',strNetwork,strNetwork);
savedata1(path,[Thirty_Run_maxNMI;0 0;mean(Thirty_Run_maxNMI(:,1)) mean(Thirty_Run_maxNMI(:,2))]); 
